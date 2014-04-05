/* fisherCMB.cc
 * Contains functions for calculating fisher matrix for a CMB experiment
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "fisherCMB.h"
#include "fisher.h"
#include "dcosmology.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include "dnumrecipes.h"
#include "spline.h"

using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 FisherCMB::FisherCMB()
{
	//cout<<"FisherCMB constructor called"<<endl;
	nchannel=0;
	fisher_flag=0;
        nparamFC=nparamCORE+nparamCMBX;
       
}

 FisherCMB::~FisherCMB()
{
	//cout<<"FisherCMB destructor called"<<endl;

	if(nchannel>0){
		free_dvector(beam,1,nchannel);
		free_dvector(sigmat,1,nchannel);
		free_dvector(sigmap,1,nchannel);
	}	
}

//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////

//initialise CMB experiment parameters
void FisherCMB::initExperiment(string name1, double fsky1, int xmax1, int lmax1, double nu1, double beam1, double sigmat1, double sigmap1)
{

	exp_name=name1;
	fsky=fsky1;
	xmax=xmax1;
	lmax=lmax1;
	addChannelExperiment(nu1,beam1,sigmat1,sigmap1);

}

void FisherCMB::addChannelExperiment(double nu1, double beamin, double sigmat1, double sigmap1)
{
	double **beamsave;
        double beam1,arcmin(2.90888e-4); //conversion between arcmin and radians
	int i;

	beam1=beamin*arcmin;

	//increment channel counter
	nchannel++;

	if(nchannel>1){
		//store old beam data
		beamsave=dmatrix(1,nchannel-1,1,3);
		for(i=1;i<=nchannel-1;i++){
			beamsave[i][1]=beam[i];
			beamsave[i][2]=sigmat[i];
			beamsave[i][3]=sigmap[i];			
		}
		//free up old beam vectors
		free_dvector(beam,1,nchannel-1);
		free_dvector(sigmat,1,nchannel-1);
		free_dvector(sigmap,1,nchannel-1);

		//resize and enter data
		beam=dvector(1,nchannel);
		sigmat=dvector(1,nchannel);
		sigmap=dvector(1,nchannel);

		//enter old data
		for(i=1;i<=nchannel-1;i++){
			beam[i]=beamsave[i][1];
			sigmat[i]=beamsave[i][2];
			sigmap[i]=beamsave[i][3];			
		}

		//enter new data
		beam[nchannel]=beam1;
		sigmat[nchannel]=sigmat1;
		sigmap[nchannel]=sigmap1;		

		//clean up 
		free_dmatrix(beamsave,1,nchannel-1,1,3);
	}else{
		beam=dvector(1,nchannel);
		sigmat=dvector(1,nchannel);
		sigmap=dvector(1,nchannel);
		beam[nchannel]=beam1;
		sigmat[nchannel]=sigmat1;
		sigmap[nchannel]=sigmap1;
	}

}

// change info for a channel
int FisherCMB::replaceChannel(int nreplace, double nu1, double beam1, double sigmat1, double sigmap1)
{
	if(nreplace>nchannel){
		cout<<"optimist!  Not enough channels"<<endl;
		return 0;
	}

	beam[nreplace]=beam1;
	sigmat[nreplace]=sigmat1;
	sigmap[nreplace]=sigmap1;

	return 1;
}

// change polarization details
void FisherCMB::changePolarization(int xmax1)
{
	xmax=xmax1;

	 if(xmax==1)  cout <<"CMB Errors: Temperature only"<<endl;
 	 if(xmax==2)  cout <<"CMB Errors: T+TE"<<endl;
 	 if(xmax==3)  cout <<"CMB Errors: T+TE+EE"<<endl;
  	 if(xmax==4)  cout <<"CMB Errors: T+TE+EE+BB"<<endl;

}

///////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////
void FisherCMB::setLMax(int lmax1)
{
	lmax=lmax1;
}



////////////////////////////////////////////////////////////////////////
// CMB Fisher matrix code
///////////////////////////////////////////////////////////////////////

//Evaluate the CMB fisher matrix 
void FisherCMB::fisherMatrixCMB()
{
  double fM;
  string file_dCli;
  string file_dClj;

  string namet, base, tag;
  int i,j,l,m;
  double luse;
  ifstream fini,finj;
  ifstream fin_ref;
  ofstream fout;
  double *dCli, *dClj,*Cl;
  string name=dirbase;
  string file_fid;
  CosmParam *c; 

  int nmat;
  double **fisher, **ifisher;

  c=&fiducialFC;

 file_fid=c->cl_file;

  nparamFC=param_tags.size()-1;
  string file_deriv[nparamFC];
//for(i=0;i<=nparamFC;i++) cout<<param_tags[i]<<endl;
  nmat=nparamFC;  
  fisher=dmatrix(1,nmat,1,nmat);
  ifisher=dmatrix(1,nmat,1,nmat);

  dCli=dvector(1,4);
  dClj=dvector(1,4);
  Cl=dvector(1,4);

  //specify derivative files
  tag="_cl_deriv.dat";
  for(i=0;i<=nparamFC-1;i++) file_deriv[i]=name+param_tags[i+1]+tag;

  //everything is done using files to create the fisher matrix
  //I'll exploit this to automate the procedure

  //parameter loops
  for(i=0;i<nmat;i++){
    for(j=0;j<=i;j++){
      fM=0.0;
      file_dCli=file_deriv[i];
      file_dClj=file_deriv[j];
      fini.open(file_dCli.c_str());
      finj.open(file_dClj.c_str());
      fin_ref.open(file_fid.c_str());
      //sum over l
      for(l=2;l<=lmax;l++){
	//use Cl in (\mu K)^2
	//reorder to better suit experimental possibilities
	//x=(T,E,B,C) -> x=(T,C,E,B)
	fini>>luse>>dCli[1]>>dCli[3]>>dCli[4]>>dCli[2];
	finj>>luse>>dClj[1]>>dClj[3]>>dClj[4]>>dClj[2];
	fin_ref>>luse>>Cl[1]>>Cl[3]>>Cl[4]>>Cl[2];
	//recover the Cl from l(l+1)Cl/2*PI
	for(m=1;m<=4;m++){
	  Cl[m]/=double(luse*(luse+1))/(2.0*PI);
	  dCli[m]/=double(luse*(luse+1))/(2.0*PI);
	  dClj[m]/=double(luse*(luse+1))/(2.0*PI);
	}
	fM+=formFisherCMB(l,Cl,dCli,dClj);
	//cout <<i<<"\t"<<j<<"\t"<<luse<<"\t"<<fM<<endl;
      }
      fisher[i+1][j+1]=fM;
      if(i!=j) fisher[j+1][i+1]=fM;  //fill in matrix by symmetry
      fini.close();
      finj.close();
      fin_ref.close();
      //      cout <<i<<"\t"<<j<<"\t"<<fM<<endl;
    }
  }
  free_dvector(dCli,1,4);
  free_dvector(dClj,1,4);
  free_dvector(Cl,1,4);

  //save information about fisher matrix and inverse into class
  setFisherMatrix(fisher,nmat,c);

  //now write results to file
  if(xmax==1) tag="_T.dat";
  if(xmax==2) tag="_TC.dat";
  if(xmax==3) tag="_TCE.dat";
  if(xmax==4) tag="_TCEB.dat";

  saveFisher(tag);
  
  cout<<exp_name<<endl;
  cout<<nparamFC<<endl;
  if(xmax==1)  cout <<"CMB Errors: Temperature only"<<endl;
  if(xmax==2)  cout <<"CMB Errors: T+TE"<<endl;
  if(xmax==3)  cout <<"CMB Errors: T+TE+EE"<<endl;
  if(xmax==4)  cout <<"CMB Errors: T+TE+EE+BB"<<endl;

  for(i=1;i<=nparamFC;i++){
	cout<<"delta"<<param_tags[i]<<"\t"<<param_values[i]<<"\t"<<sqrt(ifisherFC[i][i])<<"\t"<<sqrt(ifisherFC[i][i])/param_values[i]<<endl;
  }

  free_dmatrix(fisher,1,nmat,1,nmat);
  free_dmatrix(ifisher,1,nmat,1,nmat);
}


double FisherCMB::formFisherCMB(int luse,double *Cl,double *dCli,double *dClj)
{
  double fFC(0.0);
//  int xmax;
  int x,y;
  double **cov,**icov;

  cov=dmatrix(1,xmax,1,xmax);
  icov=dmatrix(1,xmax,1,xmax);

  //form covariance matrix
  for(x=1;x<=xmax;x++){
    for(y=1;y<=xmax;y++){
      cov[x][y]=covCMB(luse,x,y,Cl);
    }
  }
  //obtain inverse covariance matrix
  invertMatrix(cov,xmax,icov);

  //loop over polarisation types
  for(x=1;x<=xmax;x++){
    for(y=1;y<=xmax;y++){
      fFC+=dCli[x]*dClj[y]*icov[x][y];
    }
  }
  return fFC;
}

//covariance matrix for a CMB experiment
//x=(T,C,E,B)
double FisherCMB::covCMB(int l,int x,int y,double *Cl)
{
  double cov;
  double bl2;
  double bl2temp;
  double ll;
  double wtbl2, wpbl2;
  static double iwtbl2,iwpbl2;
  static double prefactor;
  static int lsave;
	int i;

  ll=double(l);

  //calculate weighting functions anew for each value of l
  if(lsave !=l){
    lsave=l;
	//sum contribution from all channels
	wtbl2=0.0;
	wpbl2=0.0;
	for(i=1;i<=nchannel;i++){
    		bl2temp=pow(beam[i],2.0)*ll*(ll+1.0);
    		bl2temp/=-8.0*log(2.0);
    		bl2=exp(bl2temp);
    		wtbl2+=pow(sigmat[i]*beam[i],-2.0)*bl2;
    		wpbl2+=pow(sigmap[i]*beam[i],-2.0)*bl2;
	}	
    iwtbl2=1.0/wtbl2;
    iwpbl2=1.0/wpbl2;
    prefactor=2.0/(2.0*ll+1.0)/fsky;
  }

  //evaluate covariance matrices
  //x=(T,C,E,B)  note this is different from CAMB output ordering

  cov=prefactor;
  if(x==1 && y==1) cov*=pow(Cl[1]+iwtbl2,2.0);
  if(x==1 && y==3) cov*=pow(Cl[2],2.0);
  if(x==1 && y==2) cov*=Cl[2]*(Cl[1]+iwtbl2);
  if(x==1 && y==4) return 0.0;
  if(x==3 && y==1) cov*=pow(Cl[2],2.0);  
  if(x==3 && y==3) cov*=pow(Cl[3]+iwpbl2,2.0); 
  if(x==3 && y==2) cov*=Cl[2]*(Cl[3]+iwpbl2);
  if(x==3 && y==4) return 0.0;
  if(x==2 && y==1) cov*=Cl[2]*(Cl[1]+iwtbl2);    
  if(x==2 && y==3) cov*=Cl[2]*(Cl[3]+iwpbl2);
  if(x==2 && y==2) cov*=(pow(Cl[2],2.0)+(Cl[1]+iwtbl2)*(Cl[3]+iwpbl2))/2.0;
  if(x==2 && y==4) return 0.0;
  if(x==4 && y==1) return 0.0;  
  if(x==4 && y==3) return 0.0;
  if(x==4 && y==2) return 0.0;
  if(x==4 && y==4) cov*=pow(Cl[4]+iwpbl2,2.0);

  return cov;
}


