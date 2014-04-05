

/*  reionization.cc
 *
 *  contains code for reproducing the Furlanetto, Zaldarriaga, Hernquist
 *  bubble model for expanding HII regions.
 *
 *  contains other code useful for reionization
 *
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include "astrophysics.h"
#include "dcosmology.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "dnumrecipes.h"
#include "ionization.h"
#include "haloDensity.h"
#include "reionization.h"
#include <sstream>


using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

Ionization::Ionization(Cosmology *c1, Astrophysics *a1)
{
  double mb(1.0),dsdmB;
  // cout <<"Reionization constructor has been called." <<endl;
  c=c1;
  a=a1;
  zion_cut=ZION_START;
  ZREION=-1.0;
  clumping_param=1.0;
  globalHistoryFlag=0;
  sigm(mb,dsdmB,c,2);

  cmb_flag=0;
  ndot_flag=0;
  lls_flag=1;
  nocut_flag=0;
  nocmb_flag=0;

  nthread=0;
  filetag=".dat";
  pdfdir="";
  gamma_flag=1;
  tocm_flag=0;
}

Ionization::~Ionization()
{
  // cout <<"Reionization destructor has been called." <<endl;

}

/////////////////////////////////////////////////////////////////////
// Member functions
/////////////////////////////////////////////////////////////////////
void Ionization::setFlags(int ndot, int cmb, int lls, int nt)
{
  stringstream ss;
   ndot_flag=ndot;
   cmb_flag=cmb;
   lls_flag=lls;
   nthread=nt;
   if(nthread>0){
     //sprintf(file,"_%i.dat",nthread);
     ss<<"_"<<nthread<<".dat";
     ss>>filetag;
     cout<<filetag<<endl;
   }


   if(cmb_flag==0) cout<<"using wmap5"<<endl;
   if(cmb_flag==1) cout<<"using wmap3"<<endl;
   if(cmb_flag==2) cout<<"using planck"<<endl;

   //if(param_flag==2) cout<<"using twostep"<<endl;
   //if(param_flag==3) cout<<"using ndot"<<endl;

   if(ndot_flag==0) cout<<"using bolton ndot"<<endl;
   if(ndot_flag==1) cout<<"using fg ndot"<<endl;
   if(ndot_flag==2) cout<<"ignoring Lya forest"<<endl;

   if(lls_flag==0) cout<<"assuming no lls correction"<<endl;
}

double Ionization::getZeta(double z)
{
  vector<double> empty;
  return setZetaParam(z,empty,0);
}

double Ionization::getZCut()
{
  return zion_cut;
}

double Ionization::getClumpingParam()
{
	return clumping_param;
}

double Ionization::setZetaParam(double z, vector<double> zparamin, int iflag)
{
   static int init_flag(0);
	double z0(4.0);
	double zeta(1.0);
	double exp_zeta(0.0);
	double NMax(10.0);
	int i;
	static vector<double> zparam;
	static double gamma4;
	double tol(1.0e-4);
	  double dzdt,Nion,NionF;
	  double Nion0(1.0e51/MPC/MPC/MPC);  //normalization

	//will need to recalculate history
	//globalHistoryFlag=0;

	if(iflag==1){
           zparam.clear();
           for(i=0;i<=zparamin.size()-1;i++) zparam.push_back(zparamin[i]);
           globalHistoryFlag=0; //need new history for new param
           if(zparam[0]==3)  zion_cut=zparam[2];
           if(zparam[0]==4){
	     zion_cut=zparam[2];
	     setDummyGamma(gamma4,zparam[2]-z0,zparam[3], zparam[4],1);
	     gamma4=getGamma4();
	     zion_cut=fmin(checkPolynomial(z0,zparam[3],zparam[4],gamma4),zion_cut);
	     if(checkNionMax(zparam[1],zparam[3],zparam[4],gamma4,zion_cut)>NMax) zion_cut=1.0;
	     //zion_cut=fmax(zion_cut,7.0);
	     //cout<<zion_cut<<"\t"<<zparam[2]<<endl;
	     //cout<<zparam[1]<<"\t"<<zparam[2]<<"\t"<<zparam[3]<<"\t"<<zparam[4]<<"\t"<<gamma4<<"\t"<<polynomial(zion_cut-z0,zparam[3],zparam[4],gamma4)<<"\t"<<zion_cut<<"\t"<<zparam[2]<<endl;
	   }
	   
           init_flag=1;
           return 0.0;
	}


        if(init_flag==0) zeta=40.0;

	//constant zeta parametrization
	if(zparam[0]==0) return zparam[1];

	if(zparam[0]==1){
	//exponential parameterization
	zeta=zparam[1];
	if(zparam.size()>2) exp_zeta+=zparam[2]*(z-z0);
	if(zparam.size()>3) exp_zeta+=zparam[3]*pow(z-z0,2.0);
	if(zparam.size()>4) exp_zeta+=zparam[4]*pow(z-z0,3.0);
	zeta*=exp(exp_zeta);
	

	}else if(zparam[0]==2){
	//tanh parametrization
	zeta=zparam[1];
	if(zparam.size()>4){
	   zeta+=(zparam[2]-zparam[1])*(tanh((z-zparam[3])/zparam[4])+1.0)/2.0;
	}
	
	}else if(zparam[0]==3){
	  //polonomial for Nion
	  z0=4.0;  //match to lowest data point
	  dzdt=(1.0+z)*c->hubbleZ(z)*UNH*c->getH();
	  NionF=c->nh(0.0)*a->dfcdzFast(z)*dzdt;  //correct for zeta param
	  Nion=Nion0*zparam[1];
	  if(zparam.size()>2) zion_cut=zparam[2];
	  if(zparam.size()>3) exp_zeta+=zparam[3]*(z-z0);
	  if(zparam.size()>4) exp_zeta+=zparam[4]*pow(z-z0,2.0);
	  if(zparam.size()>5) exp_zeta+=zparam[5]*pow(z-z0,3.0);
	  Nion*=1.0+exp_zeta;
	  if(Nion<0.0){
	    cout<<"neg: "<<zparam[1]<<"\t"<<zparam[2]<<"\t"<<zparam[3]<<"\t"<<zparam[4]<<"\t"<<z<<"\t"<<Nion<<endl;
	    Nion=1.0e-10;
	  }
	  zeta=Nion/NionF;
	}else if(zparam[0]==4){
	  //polynomial that crosses Nion=0 at xcut

	  z0=4.0;  //match to lowest data point
	  dzdt=(1.0+z)*c->hubbleZ(z)*UNH*c->getH();
	  NionF=c->nh(0.0)*a->dfcdzFast(z)*dzdt;  //correct for zeta param
	  Nion=Nion0*zparam[1]*polynomial(z-z0,zparam[3],zparam[4],gamma4);

	  if(Nion<0.0){
	    //cout<<"neg: "<<zparam[1]<<"\t"<<zparam[2]<<"\t"<<zparam[3]<<"\t"<<zparam[4]<<"\t"<<gamma4<<"\t"<<z<<"\t"<<Nion<<endl;
	    Nion=1.0e-4*Nion0;
	  }
	  zeta=Nion/NionF;
	  
	}else if(zparam[0]==5){
	//polynomial in zeta parametrization
	  zeta=zparam[1]*polynomial(z-z0,zparam[2],zparam[3],zparam[4]);       
	  if(zeta<0.0) zeta=0.1;
	}


	//cout<<z<<"\t"<<zeta<<endl;

	return zeta;
}

double polynomial(double x, double a, double b, double c)
{
  return 1.0+a*x+b*x*x+c*x*x*x;
}

double checkPolynomial(double z0, double a, double b, double c)
{
  double z;
  double zstep(0.001);
  double poly(1.0);

  z=z0;
  zstep=0.1;
  while(poly>0.0 && z<=40.0){
    poly=polynomial(z-z0,a,b,c);
    z+=zstep;
  }
  z-=2.0*zstep;
  //cout<<"polynomial0= "<<z<<"\t"<<polynomial(z-z0,a,b,c)<<"\t"<<polynomial(z+zstep-z0,a,b,c)<<endl;

  zstep=0.001;
  poly=polynomial(z-z0,a,b,c);
  while(poly>0.0 && z<=40.0){
    poly=polynomial(z-z0,a,b,c);
    z+=zstep;
  }
  z-=2.0*zstep;
  //cout<<"polynomial1= "<<z<<"\t"<<polynomial(z-z0,a,b,c)<<"\t"<<polynomial(z+zstep-z0,a,b,c)<<endl;

  zstep=0.00001;
  poly=polynomial(z-z0,a,b,c);
  while(poly>0.0 && z<=40.0){
    poly=polynomial(z-z0,a,b,c);
    z+=zstep;
  }
  z-=2.0*zstep;

  if(polynomial(z-z0,a,b,c)<0.0){
    cout<<"Error in checkPolynomial: still negative at z="<<z<<endl;
    cout<<"param: "<<a<<"\t"<<b<<"\t"<<c<<endl;
  }
  //cout<<"polynomial3= "<<z<<"\t"<<polynomial(z-z0,a,b,c)<<"\t"<<polynomial(z+zstep-z0,a,b,c)<<endl;

  return z;
}

double Ionization::getGamma4()
{
  double tol(1.0e-4);
  double cmin(-0.2), cmax(0.2);
  int icount(0);
  double gmin, gmax;
  double gamma4;

  icount=0;
  gmin=dummyGamma(cmin);
  gmax=dummyGamma(cmax);
  while(fabs(gmin/fabs(gmin)-gmax/fabs(gmax))<1.0e-4){
    //cout<<"improving gamma4"<<endl;
    //cout<<gmin<<"\t"<<gmax<<endl;
    if(icount%2==0){
      cmin-=0.05;
    }else{
      cmax+=0.05;
    }
    gmin=dummyGamma(cmin);
    gmax=dummyGamma(cmax);
    //cout<<cmin<<"\t"<<cmax<<"\t"<<gmin<<"\t"<<gmax<<"\t"<<icount<<endl;
    icount++;
  }
  gamma4=zriddrSimp(dummyGamma,cmin,cmax,tol);
  //cout<<gamma4<<endl;
  return gamma4;
}

double dummyGamma(double gamma)
{
  return setDummyGamma(gamma,0.0,0.0,0.0,0);
}

double setDummyGamma(double gamma, double zcut1, double a1, double b1, int iflag)
{
  static double a, b, zcut;

  if(iflag==1){
    a=a1;
    b=b1;
    zcut=zcut1;
    return 0.0;
  }
  return polynomial(zcut,a,b,gamma);
}

//find maximum value of Nion within integration range
double Ionization::checkNionMax(double Nion, double a, double b, double c, double zcut)
{
  double nuse, nmax(-1.0e30);
  double z;
  double z0(4.0);

  z=z0;
  while(z<=zcut){
    nuse=Nion*polynomial(z-z0,a,b,c);
    nmax=fmax(nuse,nmax);
    z+=0.1;
  }

  return nmax;
}

double Ionization::getNion(double z)
{
   double dzdt,nion;
   double step;
   double dz(0.05);

   //odeint complains if switch on background abruptly
   step=1.0-(tanh((z-zion_cut)/dz)+1.0)/2.0;

   dzdt=(1.0+z)*c->hubbleZ(z)*UNH*c->getH();
   nion= c->nh(0.0)*getZeta(z)*a->dfcdzFast(z)*dzdt;
   //cout<<z<<"\t"<<getZeta(z)<<"\t"<<nion<<endl;
   return nion*step;
}

//set pdf for ndot likelihoods
void Ionization::setPDF(string file)
{
  pdfdir=file;
}

void Ionization::setGammaFlag(int gamma_flag1)
{
  gamma_flag=gamma_flag1;
}

void Ionization::setTocmFlag(int tocm1)
{
  tocm_flag=tocm1;
}

void Ionization::setNocutFlag(int nocut1)
{
  nocut_flag=nocut1;
}

void Ionization::setNocmbFlag(int nocmb1)
{
  nocmb_flag=nocmb1;
}

/////////////////////////////////////////////////////////////////////
// Integrate ionized fraction evolution
/////////////////////////////////////////////////////////////////////

// function to calculate derivatives for globalHistory
void setDerivsHistoryIon(double z, double y[],double dy[],Cosmology *c1,Astrophysics *a1, Ionization *Ion1, int iflag)
{
  static Cosmology *c;
  static Astrophysics *a;
  static Ionization *Ion;
  double dx;
  double dzdt,zz;
  double xi;
  //double zeta;  // ionisation efficiency parameter
  double ion,recomb;
  double alphaA(4.2e-13);  //case-A recombination: cm^3 s^-1 (T=10^4K)
  double C(1.0);  //Clumping factor


  if(iflag==1){
    c=c1;
    a=a1;
    Ion=Ion1;
  }else{

    xi=y[1];
    if(xi>1.0) xi=1.0;
    if(xi<0.0) xi=0.0;

    if(clumping_flag==1) C=Ion->getClumping(z,xi);

    //////////////////////////////
    //Handle ionization evolution, x_i
    zz = z+1.0;
    dzdt=zz*c->hubbleZ(z)*UNH*c->getH();

    ion=Ion->getNion(z)/c->nh(0.0)/dzdt;

    recomb=alphaA*C*xi*c->nh(z);  
    recomb/=dzdt;
    
    dx=ion-recomb;
    
    dy[1]= -dx;
    if(fabs(xi-1.0)<5.0e-4){
    //if(xi-1.0> -1.0e-6){
      dy[1]=0.0; //cheat in case that xfree>1 
      //y[1]=1.0;
    }
    //cout<<"step: "<<z<<"\t"<<xi<<"\t"<<C<<"\t"<<dx<<"\t"<<y[1]<<"\t"<<dy[1]<<endl;
   }
}

//dummy function to get derivative for odeint for globalHistory
void derivsHistoryIon(double z, double y[], double dy[])
{
  setDerivsHistoryIon(z,y,dy,NULL,NULL,NULL,0);
}

////////////////////////////////////////////////////////////////
// Spline for getting xi quickly
////////////////////////////////////////////////////////////////

// Integrate up the history of T_IGM and x_e
// Output full history in "Thistory_save.dat"
// Return values of T and x_e at zin
//
//Units: T  K
void Ionization::globalHistory(void)
{
  double z,zstart(zion_cut);
  double zstep(0.1), zmin(ZION_END);
  double *y;
  double *dy;
  //double globalxi;
  int nvar(1);
  int nok,nbad;
  double hstep(0.1);
  double hmin(0.0001);
  double eps(1e-4);

  ofstream fout;
  string file;

  double *zs, *xis;
  double datastore[500][2];
  int i,ndata(1);

  zstart=zion_cut;

  y=dvector(1,nvar);
  dy=dvector(1,nvar);

  //cout<<"starting"<<endl;
  
  if(clumping_flag==1) setClumping();

  setDerivsHistoryIon(zmin,y,dy,c,a,this,1);

  //cout<<"calculating global history"<<endl;

  //call RECFAST and extract T_IGM and x_e history before reionization

  y[1]=1.0e-5;  //initially there are no bubbles
  //if(zstart<7.0) y[1]=0.01;  //ugly hack don't like at all

  //store initial conditions then loop over decreasing redshift
  datastore[ndata][0]=zstart;
  datastore[ndata][1]=y[1];
  while(zstart>zmin && zstart>0.0){
    ndata++;
    z=zstart-zstep;
    //cout<<z<<"\t"<<y[1]<<endl;
    if(y[1]<1.0-5.0e-4){
    odeint(y,nvar,zstart,z,eps,hstep,hmin,&nok,&nbad,derivsHistoryIon,rkqs);
    if(y[1]>1.0) y[1]=1.0;
    }else{
      y[1]=1.0;
    }

    datastore[ndata][0]=z;
    datastore[ndata][1]=y[1];
    //cout<<z<<"\t"<<y[1]<<endl;
    zstart=z;
  }

  zs=dvector(1,ndata);
  xis=dvector(1,ndata);

  // output calculated information in ascending redshift order
  //file="Thistory_save"+filetag;
  //cout<<file<<endl;
  //fout.open(file.c_str());
  for(i=ndata;i>0;i--){
    //fout<<datastore[i][0]<<"\t"<<datastore[i][1]<<endl;

    zs[ndata-i+1]=datastore[i][0];
    xis[ndata-i+1]=datastore[i][1];
  }
  //fout.close();

  //set splines for future use
  globalXI.setSplineSP(ndata,zs,xis);

  //tidy memory
  free_dvector(y,1,nvar);  
  free_dvector(dy,1,nvar);
  free_dvector(zs,1,ndata);
  free_dvector(xis,1,ndata);

  //cout<<"global history calculated successfully"<<endl;
  globalHistoryFlag=1;
  ZREION=findZReion();
  //cout<<"reionization at z="<<ZREION<<endl;
}



double Ionization::getXI(double z)
{
  if(zion_cut<6.0) return 0.0;
  if(globalHistoryFlag==0) globalHistory();

  if(z>zion_cut) return 0.0;
  if(z<ZION_END) return 1.0;
  if(z<ZREION) return 1.0;

  return globalXI.returnValue(z);
}

// From global history find out when the Universe became ionized
//
double Ionization::findZReion()
{
  ifstream fin;
  string file;
  double z, xi;
  //double tk, xe;
  double zmin;
  double zmark(-1.0);

  //use file
  /*
  file="Thistory_save"+filetag;
  fin.open(file.c_str());
  while(!fin.eof()){
    fin>>z>>xi;
    //    cout<<z<<"\t"<<tk<<"\t"<<xi<<"\t"<<xe<<endl;
    if(xi>1.0-1.0e-4) zmark=z;
  }
  */

  //use spline
  z=zion_cut;
  zmin=globalXI.getXMin();
  while(z>zmin){
    xi=globalXI.returnValue(z);
    //    cout<<z<<"\t"<<tk<<"\t"<<xi<<"\t"<<xe<<endl;
    if(xi>1.0-1.0e-4) zmark=z;
    z-=0.1;
  }

  if(zmark<0.0)  zmark=ZION_END;

  return zmark;
}


////////////////////////////////////////////////////////////////
// Calculate tau for a specified ionization history
////////////////////////////////////////////////////////////////

//get optical depth to CMB
double Ionization::getTauCMB()
{
  double tau;
  double zmin(0.0),zmax(ZION_START);
  double tol(1.0e-4);
  double zreion(ZION_END);
  double zHe(3.0);

  zmax=zion_cut;

  setTauIntIon(zmin,zreion,c,this,0);
  tau=qromb(getTauIntIon,zHe+1.0e-4,zmax,tol);
  tau+=qromb(getTauIntIon,zmin,zHe-1.0e-4,tol);

  return tau;
}

//integrand for tau integral
double setTauIntIon(double z, double zri1,Cosmology *c1, Ionization *Ion1,int iflag)
{
  double dtaudz;
  static Cosmology *c;
  static Ionization *Ion;
  static double zri;
  double zHe(3.0);
  double xi;

  if(iflag==0){
    zri=zri1;
    c=c1;
    Ion=Ion1;
    return 0.0;
  }

  dtaudz=c->nh(z)*SPEEDOFLIGHT_CGS*SIGMATHOMSON;
  dtaudz/=(1.0+z)*c->hubbleZ(z)*UNH*c->getH();
  if(z>zri){
    xi=Ion->getXI(z)*(1.0+FHE);   //account for HeI ionization
    dtaudz*=xi; 
  }
  if(z<zHe) dtaudz*=1.0+2.0*FHE;  //account for HeII ionization

  return dtaudz;
} 

//dummy function for tau integral
double getTauIntIon(double z)
{
  return setTauIntIon(z,0.0,NULL,NULL,1);
}

////////////////////////////////////////////////////////////////
// Likelihood calculation
////////////////////////////////////////////////////////////////

double Ionization::likelihoodZeta(vector<double> param)
{
   //int cmb_flag(0);   //0 WMAP5  1 WMAP3  2 PLANCK
   //int ndot_flag(0);  //0 BOLTON 1 FG
   double z;
   double gamma3, gamma4, gamma5, gamma6;
   double tau;
   double fudge(1.0);
   int nogamma_flag(0);
   
   double tauW(0.087), sigW(0.017);
   double gammaObs3(0.86), sigGammaObs3(0.3);
   double gammaObs4(0.97), sigGammaObs4(0.4);
   double gammaObs5(0.52), sigGammaObs5(0.3);
   double gammaObs6(0.19), sigGammaObs6(0.19);
     double z_tocm(9.5);
     double xi_tocm(0.5);
     double sig_tocm(0.05);
   static int initrun(1);
   string file;

   double like(1.0);
   Lyman lyf;
   lyf.initLyman(c,a);
   lyf.setLLSFlag(lls_flag);
 
   setZetaParam(0.0,param,1);
   if(zion_cut<6.0) return 0.0;  
   //don't believe any model like this since it won't reionize early enough
   //these models cause problems in global history, so catch them here.

   //cout<<getXI(6.0)<<endl;
   if(nocut_flag==0){
     if(getXI(6.0)<0.99) return 0.0;
   }

   if(cmb_flag==1){
      //WMAP3
      tauW=0.09;
      sigW=0.03;
   }else if(cmb_flag==2){
      //PLANCK
      tauW=0.09;
      //sigW=0.01;  //includes likely effect of foreground modeling
      sigW=0.005; //no foregrounds
   }

   if(ndot_flag==1){
      //fg constraints
      gammaObs4=0.55;
      sigGammaObs4=0.05;
   }

   
   //CMB optical depth
   if(nocmb_flag==0){
     tau=getTauCMB();
     like=gaussianProb(tau,tauW,sigW);
   }

   //cout<<tau<<endl;
   if(ndot_flag==2) nogamma_flag=1;

   if(nogamma_flag==0){

   if(gamma_flag==1){
     // calculate ionization rate
     z=4.0;
     gamma4=lyf.getGammaFromNdot(z,getNion(z)*MPC*MPC*MPC);
     z=5.0;
     gamma5=lyf.getGammaFromNdot(z,getNion(z)*MPC*MPC*MPC);
     z=6.0;
     gamma6=lyf.getGammaFromNdot(z,getNion(z)*MPC*MPC*MPC);
  
     like*=gaussianProb(gamma4,gammaObs4,sigGammaObs4*fudge);

     if(ndot_flag==0){
       like*=gaussianProb(gamma5,gammaObs5,sigGammaObs5*fudge);
       //like*=gaussianProb(gamma6,gammaObs6,sigGammaObs6*fudge);
     }

   }else{

     z=4.0;
     file=pdfdir+"zpdf_astro_z4.dat";
     like*=nongaussianProb(getNion(z)*MPC*MPC*MPC/1.0e51,&ndot4sp,file,initrun);

     z=5.0;  file=pdfdir+"zpdf_astro_z5.dat";
       like*=nongaussianProb(getNion(z)*MPC*MPC*MPC/1.0e51,&ndot5sp,file,initrun);


     if(ndot_flag==0){
       z=6.0;  file=pdfdir+"zpdf_astro_z6.dat";
       like*=nongaussianProb(getNion(z)*MPC*MPC*MPC/1.0e51,&ndot6sp,file,initrun);
     }
     initrun=0;
   }
   }//end nogamma if

   //21 cm
   if(tocm_flag>0){
     z_tocm=9.0;
     xi_tocm=0.5;
     sig_tocm=0.05;
     like*=gaussianProb(getXI(z_tocm),xi_tocm,sig_tocm);
   }
   if(tocm_flag>1){
     z_tocm=7.0;
     xi_tocm=0.9;
     sig_tocm=1.0;
     like*=uniformProb(getXI(z_tocm),xi_tocm,sig_tocm);     
   }


   if(like<1.0e-40) like=0.0;  //values of 1e-300 cause stream errors

   return like;
}

//Normalised Gaussian probability function
double Ionization::gaussianProb(double x, double x0, double sig)
{
  double prob;
  
  prob=exp(-pow(x-x0,2.0)/2.0/sig/sig);
  prob/=sqrt(2.0*PI)*sig;

  return prob;
}

double Ionization::uniformProb(double x, double xmin, double xmax)
{
  double width;
  width=(xmax-xmin);

  if(x<xmin || x>xmax) return 0.0;
  
  return 1.0/width;
}

//obtain probability from spline - retains non-gaussian likelihood
double Ionization::nongaussianProb(double x, Spline *spline, string file, int iflag)
{
  if(iflag==1){
    spline->loadFileSpline(file,3,1,3);
  }

  //cout<<x<<"\t"<<spline->getXMin()<<"\t"<<spline->getXMax()<<endl;
  if(x<spline->getXMin()) return 0.0;
  if(x>spline->getXMax()) return 0.0;
  return spline->returnValue(x);
}

////////////////////////////////////////////////////////////////
// Two dimensional spline for clumping to speed up calculation
////////////////////////////////////////////////////////////////

// Since we're most interested in the behaviour of the clumping when
// x_H is close to zero the spline works better if we spline in x_H rather
// than x_i.  This allows us to get arbitrarily close to x_H=0.  The spline
// is less well behaved at x_H=1, but since the clumping is negligible there
// this should not affect my calculations significantly.
//
void Ionization::setClumping()
{
	int nxi(160);
	int nz(40);
	//nz=60;
	int i,j;
	double dxi, dz;
	double ximin(log(1.0e-5));
	double ximax(log(1.0-1.0e-3));
	//ximin=0.001;
	//ximax=1.0-1.0e-3;
	double zmax(ZION_START);
        zmax=29.0;
	//zmax=ZION_START;
	double zmin(ZION_END);
	double *zl, *xil, **cl;	

	ifstream fin;
	char *file="clumping_spline.dat";

	fin.open(file);

	if(fin.fail()){
	  Reionization rz(c,a);
	   ofstream fout;
	   
           fout.open(file);

           dxi=(ximax-ximin)/(double)(nxi-1);
           dz=(zmax-zmin)/(double)(nz-1);

           zl=dvector(1,nz);
           xil=dvector(1,nxi);
           cl=dmatrix(1,nz,1,nxi);

           cout<<"initialising clumping spline"<<endl;

           fout<<nz<<"\t"<<nxi<<endl;
           for(i=1;i<=nz;i++){
              zl[i]=zmin+(double)(i-1)*dz;
              a->initPVParamMEHR(zl[i]);
              for(j=1;j<=nxi;j++){
                 xil[j]=ximin+(double)(j-1)*dxi;
                 //cl[i][j]=4.0*a->clumpingMEHR(zl[i],xil[j]);
                 //cl[i][j]=4.0*a->clumpingMEHR(zl[i],1.0-exp(xil[j]));

                 cl[i][j]=rz.clumpingBub(zl[i],1.0-exp(xil[j]));
                 //cl[i][j]=rz.clumpingBub(zl[i],xil[j]);
                 cout<<zl[i]<<"\t"<<exp(xil[j])<<"\t"<<cl[i][j]<<endl;
                 //cout<<zl[i]<<"\t"<<xil[j]<<"\t"<<cl[i][j]<<endl;
                 fout<<zl[i]<<"\t"<<xil[j]<<"\t"<<cl[i][j]<<endl;
              }
           }
           fout.close();

           clumpingSP.setSplineSP(nz,nxi,zl,xil,cl);

           free_dvector(zl,1,nz);
           free_dvector(xil,1,nxi);
           free_dmatrix(cl,1,nz,1,nxi);

	}else{
           fin.close();
           clumpingSP.loadFileSpline(file);
	}

	//cout<<"clumping spline initialised"<<endl;

//	cout<<clumpingSP.getX1Min()<<"\t"<<clumpingSP.getX1Max()<<"\t"<<clumpingSP.getX2Min()<<"\t"<<clumpingSP.getX2Max()<<"\t"<<clumpingSP.getM()<<"\t"<<clumpingSP.getN()<<endl;

}

double Ionization::getClumping(double z, double xi)
{
   double clump;

   xi=1.0-xi;

   if(z>clumpingSP.getX1Max()) return 0.0;
   if(z<clumpingSP.getX1Min()){
      if(xi<exp(clumpingSP.getX2Max()) && xi>exp(clumpingSP.getX2Min()) ){
         clump=clumpingSP.returnValue(clumpingSP.getX1Min(),log(xi));
      }else if(xi<exp(clumpingSP.getX2Min())){
         clump=clumpingSP.returnValue(clumpingSP.getX1Min(),clumpingSP.getX2Min());
      }else{
         clump=clumpingSP.returnValue(clumpingSP.getX1Min(),clumpingSP.getX2Max());
      }
      return clump; 
	}
	if(xi<exp(clumpingSP.getX2Min())) return clumpingSP.returnValue(z,clumpingSP.getX2Min());
	if(xi>exp(clumpingSP.getX2Max())){
		clump=clumpingSP.returnValue(z,clumpingSP.getX2Max());
		 return clump; 
	}

	return clumpingSP.returnValue(z,log(xi));
}

//////////////////////////////////////////////////////////////////////
// Dummy function for use in VEGAS integral
//////////////////////////////////////////////////////////////////////
// use vegas as a simple of way of incorporating importance sampling
//into my code

double vegasIntegrand(double param[], double wt)
{
  string empty;
  return setVegasIntegrand(param,NULL,0,empty,0); 
}

double setVegasIntegrand(double paramin[], Ionization *Ion1, int param_flag1, string file, int iflag)
{
  static Ionization *Ion;
  static int param_flag;
  vector<double> param;
  double like, z;
  static ofstream fout;
  static int icount;

  if(iflag==1){
    Ion=Ion1;
    param_flag=param_flag1;
    fout.open(file.c_str());
    return 0.0;
  }else if(iflag==2){
    fout.close();
  }

  
  param.clear();
  param.push_back(param_flag);
  param.push_back(paramin[1]);
  param.push_back(paramin[2]);
  param.push_back(paramin[3]);
  param.push_back(paramin[4]);
  Ion->setZetaParam(0.0,param,1);
  like=Ion->likelihoodZeta(param);

 

  //handle file output
  icount++;
  cout<<icount<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like<<endl;
  fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;
  z=6.0;
  while(z<20.2){
    fout<<"\t"<<Ion->getXI(z);
    z+=0.25;
  }
  fout<<endl; 
  //end file output

  //static double pmin(1.0e30), pmax(-1.0e30);
  //if(param[1]<pmin) pmin=param[1];
  //if(param[1]>pmax) pmax=param[1];
  //cout<<pmin<<"\t"<<pmax<<endl;;

  
  //like=1.0;

  //like=pow(param[1]-15.0,2.0);
  return like;
}

