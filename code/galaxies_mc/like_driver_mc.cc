
 /* driver.cc
 * program for mimicing the results of BL04
 * Call by 
 *   driver.x -z6 > (output filename)
 *
 */

#include <iostream>
#include <fstream>
#include "math.h"
#include "astrophysics.h"
#include "dnumrecipes.h"
#include "dcosmology.h"
#include "twentyonecm.h"
#include "reionization.h"
#include "haloDensity.h"
#include "spline.h"
#include "radio.h"
#include "observation.h"
#include "fisher.h"
#include "fisherGAL.h"
#include "neutrinos.h"
#include "lymanforest.h"
#include "ionization.h"
#include "spline2D.h"

using namespace std;

/******************************************************************/


int main(int argc, char *argv[])
{
 
  // Cosmology parameters
  double om0(0.3);
  double lam0(0.7);
  double omb(0.046);
  double h(0.7);
  double s8(0.9);
  double n(1.0);
  double omNu(0.0);

  double zin_arg(15.0);
  int tin(2);
  int i,j;
  char *file;
  ofstream fout;
  ifstream fin;
  int xin(1);
  int lyaxray_in(0);
  int popflag_in(0);

     int ndot_flag(0), cmb_flag(0);
     int param_flag(0);
     int lls_flag(1);
     int like_calc_flag(1);
     int nthread(0);
     int gamma_flag(0);
     int tocm_flag(0);
     int nocut_flag(0);
     int nocmb_flag(0);
     int tau_flag(0);
     int nignore(6);


  // Handle arguments 
  while ((argc>1) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'z':
      zin_arg=atof(&argv[1][2]);
      break;
    case 't':
      nthread=atoi(&argv[1][2]);
      cout<<"thread no: "<<nthread<<endl;
      break;
    case 'x':
      xin=atoi(&argv[1][2]);
      break;
    case 'l':
      lls_flag=atoi(&argv[1][2]);
      break;
    case 'g':
      gamma_flag=atoi(&argv[1][2]);
      break;
    case 'p':
      param_flag=atoi(&argv[1][2]);
      break;
    case 'n':
      ndot_flag=atoi(&argv[1][2]);
      break;
    case 'c':
      cmb_flag=atoi(&argv[1][2]);
      break;
    case 'd':
      like_calc_flag=atoi(&argv[1][2]);
      break;
    case 'm':
      tocm_flag=atoi(&argv[1][2]);
      break;
    case 'k':
      nocut_flag=atoi(&argv[1][2]);
      break;
    case 'q':
      nocmb_flag=atoi(&argv[1][2]);
      tau_flag=1;
      tocm_flag=2;
      break;
    default:
      cerr << "Bad option" <<argv[1] <<"\n";
    }
    --argc;
    ++argv;
  }

     if(tau_flag==1) nignore++;  //now outputing cmb tau too
  
  //  s8=0.77;
  //omb=0.044;
  //om0=0.26;
  //lam0=0.74;
  //h=0.72;
  //n=0.95;

  //HQB cosmology
  //  s8=0.9;
  //omb=0.044;
  //om0=0.27;
  //lam0=0.73;
  //h=0.71;
  //n=1.0;

  //NB07 cosmology
  //  s8=0.826;
  //omb=0.0478;
  //om0=0.299;
  //lam0=0.701;
  //h=0.687;
  //n=1.0;

  //Oh1999
  //  s8=0.87;
  //omb=0.04;
  //om0=0.35;
  //lam0=0.65;
  //h=0.65;
  //n=0.96;

  cout<<"Enter cosmology parameters"<<endl;
  cout<<"Enter s8:"<<endl;
  cin>>s8;
  cout<<"Enter omb:"<<endl;
  cin>>omb;
  cout<<"Enter om0:"<<endl;
  cin>>om0;
  cout<<"Enter lam0:"<<endl;
  cin>>lam0;
  cout<<"Enter h:"<<endl;
  cin>>h;
  cout<<"Enter n:"<<endl;
  cin>>n;

  cout<<s8<<"\t"<<omb<<"\t"<<om0<<"\t"<<lam0<<"\t"<<h<<"\t"<<n<<endl;
  
  Cosmology c(om0,lam0,omb,h,s8,n,omNu);
  Astrophysics a(&c,popflag_in,xin,lyaxray_in,1.0);
  TwentyOneCM tocm(&c,&a);
  Reionization rz(&c,&a);
  Lyman lyf;

  int lyaxrayflag, popflag, xrayflag, sourceflag;
  double fstar, fx, nion, nlya, fesc;
  int modelflag;

  modelflag=tin;
  if(modelflag==1){
    //Low tau parameters
    fstar=0.1;
    fesc=0.075;
    nion=4000.0;
    nlya=-1.0;
    fx=1.0;
    lyaxrayflag=0;
    //  lyaxrayflag=0;
    popflag=0; //popII
  }else if(modelflag==2){
    //Mid tau parameters
    fstar=0.2;
    fesc=0.2;
    nion=4000.0;
    nlya=-1.0;
    fx=1.0;
    lyaxrayflag=0;
    //  lyaxrayflag=0;
    popflag=0; //popIII
  }else if(modelflag==3){
    //High tau parameters
    fstar=0.3;
    fesc=0.1;
    nion=30000.0;
    nlya=3030.0;
    fx=1.0;
    lyaxrayflag=0;
    //  lyaxrayflag=0;
    popflag=1; //popIII
  }

  xrayflag=-1;
  sourceflag=-1;

  cout<<"Enter astrophysics parameters"<<endl;
  cout<<"Enter fstar:"<<endl;
  cin>>fstar;
  cout<<"Enter fesc:"<<endl;
  cin>>fesc;
  cout<<"Enter nion:"<<endl;
  cin>>nion;
  cout<<"Enter fx:"<<endl;
  cin>>fx;
  cout<<"Enter nlya:"<<endl;
  cin>>nlya;
  cout<<"Enter popflag:"<<endl;
  cin>>popflag;
  cout<<"Enter xrayflag:"<<endl;
  cin>>xrayflag;
  cout<<"Enter lyaxrayflag:"<<endl;
  cin>>lyaxrayflag;
  cout<<"Enter sourceflag:"<<endl;
  cin>>sourceflag;

  a.initAstrophysics(fstar,fesc,nion,fx,nlya,popflag,xrayflag,lyaxrayflag,sourceflag,-1.0);

  lyf.initLyman(&c,&a);


////////////////////////////////////////////////////////////////////////
// Grid likelihood values - 1 parameter constant zeta
///////////////////////////////////////////////////////////////////////
/*
string files;
double zeta0,z;
double like;
vector<double> param;
Ionization Ion(&c,&a);
int k;
int n1(201);
double p1min,p1max,p1step;

files="likelihoodZeta_const_zeta.dat";
fout.open(files.c_str());

p1min=30.0;
p1max=70.0;

//n1=301;
//p1min=10.0;
//p1max=300.0;

p1step=(p1max-p1min)/(double)(n1-1);

for(i=1;i<=n1;i++){
	zeta0=p1min+p1step*(double)(i-1);
		
	param.clear();
	param.push_back(1);
	param.push_back(zeta0);
	param.push_back(0.0);
	param.push_back(0.0);
	param.push_back(0.0);

		Ion.setZetaParam(0.0,param,1);
		like=Ion.likelihoodZeta(param);
		cout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<like<<"\t"<<Ion.getXI(6.0)<<"\t"<<Ion.getXI(7.0)<<"\t"<<Ion.getXI(8.0)<<"\t"<<Ion.getXI(9.0)<<"\t"<<Ion.getXI(10.0)<<"\t"<<Ion.getXI(11.0)<<"\t"<<Ion.getXI(12.0)<<"\t"<<Ion.getXI(13.0)<<"\t"<<Ion.getXI(14.0)<<"\t"<<Ion.getXI(15.0)<<"\t"<<Ion.getXI(16.0)<<"\t"<<Ion.getXI(17.0)<<"\t"<<Ion.getXI(18.0)<<"\t"<<Ion.getXI(19.0)<<"\t"<<Ion.getXI(20.0)<<endl;

		fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;
		z=6.0;
		while(z<20.2){
			fout<<"\t"<<Ion.getXI(z);
			z+=0.25;
		}
		fout<<endl;

}
fout.close();

*/
////////////////////////////////////////////////////////////////////////
// Grid likelihood values - step function
///////////////////////////////////////////////////////////////////////
/*
string files;
double zeta0, zeta1,z0,deltaz;
double like,z;
double p1min, p1max, p1step;
double p2min, p2max, p2step;

vector<double> param;
Ionization Ion(&c,&a);
int k;
int n1(151);
int n2(151);

z0=13.0;
deltaz=2.0;

p1min=30.0;
p1max=75.0;
p2min=0.0;
p2max=300.0;

p1step=(p1max-p1min)/(n1-1);
p2step=(p2max-p2min)/(n2-1);

files="likelihoodZeta_twostep.dat";
fout.open(files.c_str());

for(i=1;i<=n1;i++){
	zeta0=p1min+p1step*(double)(i-1);
	for(j=1;j<=n2;j++){
		zeta1=p2min+p2step*(double)(j-1);

		param.clear();
		param.push_back(2.0);
		param.push_back(zeta0);
		param.push_back(zeta1);
		param.push_back(z0);
		param.push_back(deltaz);
		Ion.setZetaParam(0.0,param,1);
		like=Ion.likelihoodZeta(param);
		cout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<like<<endl;
		fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;
		z=6.0;
		while(z<20.2){
			fout<<"\t"<<Ion.getXI(z);
			z+=0.25;
		}
		fout<<endl;


	}
}
fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Grid likelihood values - polynomial in Nion
///////////////////////////////////////////////////////////////////////
/*
string files;
double zeta0, zeta1,z0,deltaz;
double like, z;
vector<double> param;
Ionization Ion(&c,&a);
int k;
double clump;
double p1min,p1max,p1step;
double p2min, p2max, p2step;
int ndot_flag(0), cmb_flag(0);
int param_flag;

param_flag=3;  //Nion polynomial
ndot_flag=0;
cmb_flag=0;

int n1(151);
int n2(151);
//n1=20;
//n2=20;
//n1=50;
//n2=50;

//ndot

//bolton
p1min=0.75;
p1max=1.55;

if(ndot_flag==1){
//fg
//p1min=0.83;
//p1max=1.15;
}

//cmb

p2min=6.5;
p2max=21.0;

//wmap3
if(cmb_flag==1){
p2min=6.5;
p2max=30.0;
}

p1step=(p1max-p1min)/(double)(n1-1);
p2step=(p2max-p2min)/(double)(n2-1);

deltaz=2.0;

Ion.setFlags(ndot_flag,cmb_flag);

files="likelihoodZeta_ndot.dat";
fout.open(files.c_str());

for(i=1;i<=n1;i++){
	zeta0=p1min+(double)(i-1)*p1step;
	for(j=1;j<=n2;j++){
		zeta1=p2min+(double)(j-1)*p2step;

 		param.clear();
 		param.push_back(param_flag);
 		param.push_back(zeta0);
 		param.push_back(zeta1);
 		param.push_back(0.0);
 		param.push_back(0.0);
		Ion.setZetaParam(0.0,param,1);
		like=Ion.likelihoodZeta(param);
 		cout<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like<<endl;
 		fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;
 		z=6.0;
 		while(z<20.2){
 			fout<<"\t"<<Ion.getXI(z);
 			z+=0.25;
 		}
 		fout<<endl; 

	} 
} 
fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Grid likelihood values - general case
///////////////////////////////////////////////////////////////////////
     string data_dir;
     string file_like;
     string files;
     double zeta0, zeta1,z0,deltaz;
     double zeta2, zeta3;
     double like, z;
     vector<double> param;
     Ionization Ion(&c,&a);
     //int k;
     //double clump;
     double p1, p2, p3, p4;
     double p1min(0.0),p1max(0.0),p1range;
     double p2min(0.0), p2max(0.0), p2range;
     double p3min(0.0), p3max(0.0), p3range;
     double p4min(0.0), p4max(0.0), p4range;
     int nsob,nparam(4);
     double *xsob;
     int sob_flag(1); //use sobseq rather than random points
     long idums;
     string pdfdir;

     int n1(151);
     int n2(151);
     //int npoints;
     long il, npoints;
     long ngood;
     npoints=n1*n1*n2*n2;   //number of sampling points
     npoints=n1*n2;   //~2hours to run

     npoints=500000;
     npoints=500000;
     npoints=2000000;
     //sobseq/random number initiation
     xsob=dvector(1,nparam);
     nsob=-1;
     sobseq(&nsob,xsob);  //initialisation call
     nsob=nparam;
     idums= -1;
     dran2(&idums);

     //param_flag 1: ??   2: ??  3: Nion
   //cmb_flag  0 WMAP5  1 WMAP3  2 PLANCK
   //ndot_flag 0 BOLTON 1 FG

     //param_flag=3;  //Nion polynomial
     //ndot_flag=0;
     //cmb_flag=0;

     data_dir="./";
     file_like="likelihoodZeta.dat";

     if(param_flag==2) data_dir=data_dir+"twostep/";
     if(param_flag==3 ||param_flag==4) data_dir=data_dir+"ndot/";
     if(param_flag==5) data_dir=data_dir+"polyz/";

     if(cmb_flag==1) data_dir=data_dir+"wmap3/";
     if(cmb_flag==2) data_dir=data_dir+"planck/";
     if(ndot_flag==0 && cmb_flag==0) data_dir=data_dir+"bolton/"; 
     if(ndot_flag==1) data_dir=data_dir+"fg_data/"; 
     if(ndot_flag==2) data_dir=data_dir+"nolya/"; 
     if(ndot_flag==3) data_dir=data_dir+"betterg/";

     if(lls_flag==0) data_dir=data_dir+"no_lls/";

     if(nocut_flag==1)  data_dir=data_dir+"no_cut/";
     if(nocmb_flag==1)  data_dir=data_dir+"no_cmb/";

     if(tocm_flag==1){
       data_dir=data_dir+"21cm/"; 
       cout<<"using 21 cm info"<<endl;
     }

     files="mkdir "+data_dir;
     system(files.c_str());


     if(ndot_flag==2) npoints=500000;
     /////////////////////////////////////////////////
     //set up range of parameters

     if(param_flag==2){
       z0=13.0;
       deltaz=2.0;

       p1min=30.0;
       p1max=75.0;
       p2min=0.0;
       p2max=300.0;

       p3min=6.0;
       p3max=30.0;
       p4min=0.01;
       p4max=5.0;
     }

     //////////////////////////////////////////
     if(param_flag==3){  //Nion polynomial
       //ndot

       //bolton
       p1min=0.35;
       p1max=1.55;

       if(ndot_flag==1){
	 //fg
	 p1min=0.43;
	 p1max=0.95;
       }

       //cmb

       p2min=6.5;
       p2max=30.0;

       //wmap3
       if(cmb_flag==1){
	 p2min=6.5;
	 p2max=30.0;
       }

       if(lls_flag==0){
	 p1min/=3.0;
	 p2max/=3.0;
       }

       p3min=-0.02;
       p3max=0.02;
       p4min= -0.2;
       p4max=0.2;
       p4min= -0.0005;
       p4max=0.0005;
       //p3min=-0.001;
       //p3max=0.001;
       //p4min=-0.0001;
       //p4max=0.0001;
     }

     //test code ////////////////
     if(param_flag==2){
       p1min=20.0;
       p1max=48.0;

       p2min=0.0;
       p2max=400.0;

       p3min=7.0;
       p3max=35.0;
       p4min=0.1;
       p4max=5.0;

       //p3min=13.0;
       //p3max=13.0;
       //p4min=2.0;
       //p4max=2.0;

       if(lls_flag==0){
	 p1min=0.1;
	 p1max=30.0;
	 p2min=0.0;
	 p2max=400.0;

	 p3min=6.0;
	 p3max=15.0;
	 p4min=0.1;
	 p4max=5.0;

       }

       if(ndot_flag==1 && lls_flag==1){
	 //fg
	 p1min=32.0;
	 p1min=20.0;
	 p1max=36.0;
       }

     }
     ///////////////////////////
     if(param_flag==4){  //Nion polynomial
       //ndot

       //bolton
       p1min=0.35;
       p1max=1.55;

       if(ndot_flag==1){
	 //fg
	 p1min=0.43;
	 p1max=0.95;
       }

       //cmb

       p2min=6.5;
       p2max=30.0;

       //wmap3
       if(cmb_flag==1){
	 p2min=6.5;
	 p2max=30.0;
       }

       if(lls_flag==0){
	 p1min/=3.0;
	 p2max/=3.0;
       }

       p3min=-0.1;
       p3max=0.3;
       p4min= -0.04;
       p4max=0.04;
       //p3min=-0.001;
       //p3max=0.001;
       //p4min=-0.0001;
       //p4max=0.0001;
     }
     // gamma_flag==0 cases
     //////////////////////////////////
     if(param_flag==2 && gamma_flag==0){
       p1min=10.0;
       p1max=55.0;

       p2min=0.0;

       p3min=7.0;
       p3max=35.0;
       p4min=0.1;
       p4max=5.0;

       p3max=41.0;
       p2max=505.0;
       p4max=10.1;

       p3min=4.0;
       if(ndot_flag==1){
	 //fg
	 p1min=10.0;
	 p1max=80.0;
       }

       if(ndot_flag==2) p1max=500.0;

     }

     ///////////////////////////
     if(param_flag==4 && gamma_flag==0)
       {  //Nion polynomial
       //ndot

       //bolton
       p1min=0.05;
       p1max=2.5;

       p2min=6.1;
       p2max=30.0;
       p2max=40.0;

       p3min= -2.2;
       p3max= 3.5;

       p4min= -1.2;
       p4max=2.0;

       if(ndot_flag==1){

	 p3min=-2.8;
	 p3max=2.8;

	 p4min= -1.0;
	 p4max=2.0;

	 p3min= -6.0;
	 p3max=3.5;

	 p4max=6.5;
       }

       if(ndot_flag==2) p1max=5.0;

     }

     //////////////////////////////////
     if(param_flag==5 && gamma_flag==0){
       p1min=10.0;
       p1max=65.0;

       p2min=-0.5;
       p2max=0.6;

       p3min=-0.2;
       p3max=0.4;
       p4min= -0.1;
       p4max=0.1;

       if(ndot_flag==1){
	 p1max=85.0;
	 p2max=0.8;
       }

     }

     if(tocm_flag==1){
       Ion.setTocmFlag(tocm_flag);
     }

     ////////////////////////////
     /////////////////////////////////////////////////////////
     p1range=(p1max-p1min);
     p2range=(p2max-p2min);
     p3range=(p3max-p3min);
     p4range=(p4max-p4min);

 Ion.setFlags(ndot_flag,cmb_flag,lls_flag,nthread);
 Ion.setGammaFlag(gamma_flag);
 lyf.setLLSFlag(lls_flag);
 Ion.setNocutFlag(nocut_flag);
 Ion.setNocmbFlag(nocmb_flag);

 pdfdir="pdf_data/noexternal/";
 pdfdir="pdf_data/uTuBuA/";
 pdfdir="pdf_data/modelA/";
 if(ndot_flag==1)  pdfdir="pdf_data/modelA/"; 
 if(ndot_flag==3)  pdfdir="pdf_data/modelC/";
 Ion.setPDF(pdfdir);

 if(like_calc_flag==1){

   files=data_dir+file_like;
   cout<<files<<endl;
   fout.open(files.c_str());

   //sample likelihood surface

   //for(il=1;il<=npoints;il++){

   il=0;
   while(ngood<npoints){
     il++;

     if(sob_flag==1){
       sobseq(&nsob,xsob);
       p1=p1min+xsob[1]*p1range;
       p2=p2min+xsob[2]*p2range;
       p3=p3min+xsob[3]*p3range;
       p4=p4min+xsob[4]*p4range;
     }else{
       p1=p1min+dran2(&idums)*p1range;
       p2=p2min+dran2(&idums)*p2range;
       p3=p3min+dran2(&idums)*p3range;
       p4=p4min+dran2(&idums)*p4range;
     }

     //if(il>500000){
     if(il>-1){
       //core code to get likelihood and output value
       param.clear();
       param.push_back(param_flag);
       param.push_back(p1);
       param.push_back(p2);
       param.push_back(p3);
       param.push_back(p4);
       Ion.setZetaParam(0.0,param,1);
       like=Ion.likelihoodZeta(param);
       param[2]=Ion.getZCut();  //update parameters for consistency

       if(like>0.0){ //only care about points with like>0, keeps file size down

	 ngood++;
	 if(ngood%100==0) cout<<cmb_flag<<"\t"<<ndot_flag<<"\t"<<il<<"\t"<<ngood<<"\t"<<(double)(ngood)/(double)(il)<<"\t"<<(double)(ngood)/(double)(npoints)<<endl;

	 //cout<<(double)(il)/(double)(npoints)<<"\t"<<(double)(ngood)/(double)(il)<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like<<endl;

	 fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;

	 if(tau_flag==1) fout<<"\t"<<Ion.getTauCMB(); 
	 z=6.0;
	 if(nocut_flag==1) z=4.0;
	 while(z<20.2){
	   fout<<"\t"<<Ion.getXI(z);
	   z+=0.25;
	 }
	 fout<<endl; 
       }
       //end point
     } //end if late start
   } 
   fout.close();

   // }//end if
 
 }

////////////////////////////////////////////////////////////////////////
// Process likelihood data to get neutral fraction distribution
///////////////////////////////////////////////////////////////////////

   //string data_dir;
   // string file_like;
     //string files;
     string str;
     int npoint(0);
     int ndata(0);
     int nentry;
     int nbins,nz,ibin;
     int iz;
     double pmin, pmax;
     double trash;
     //double like;
     double pdftot;
     vector<double> pv, p1v, p2v, p3v, p4v,likev, xv;
     vector<double> pdf;
     //double p1min(1e30), p1max(-1e30), p2min(1e30), p2max(-1e30);
     //char line[5];
     double likeTot(0.0);
     double xmin(0.0), xmax(0.0), dx;
     double xi,xiold,xinew;
     //double z;
     double wnorm;
     double zmin(6.0);
     double dz(0.25);
     double likeMax(0.0);
     double p0;
     double p1Best, p2Best,p3Best,p4Best;
     ofstream fout_xi;
     double evidence, pvol;
     double tau;

     if(nocut_flag==1) zmin=4.0;

    p1min=1e30;
    p1max= -1e30;
    p2min=1e30;
    p2max= -1e30;


//marginalised parameter pdfs

    files=data_dir+file_like;
    cout<<files<<endl;

    fin.open(files.c_str()); 
    while(!fin.eof()){
      fin>>pmin;
       npoint++;
       //cout<<npoint<<"\t"<<pmin<<endl;
    }
    fin.close();
    npoint--;

    cout<<npoint<<endl;

    fin.open(files.c_str()); 
    while(!fin.eof()){
       getline(fin,str);
       ndata++;
    }
    fin.close();
    ndata--;

    cout<<npoint<<"\t"<<ndata<<endl;
    //find max and min values of parameters and likelihood total
    fin.open(files.c_str()); 
    while(!fin.eof()){
       fin>>p0>>p1>>p2>>p3>>p4>>like;
       getline(fin,str);
       p1v.push_back(p1);
       p2v.push_back(p2);
       p3v.push_back(p3);
       p4v.push_back(p4);
       likev.push_back(like);

       if(p1<p1min) p1min=p1;
       if(p2<p2min) p2min=p2;
       if(p3<p3min) p3min=p3;
       if(p4<p4min) p4min=p4;
       if(p1>p1max) p1max=p1;
       if(p2>p2max) p2max=p2;
       if(p3>p3max) p3max=p3;
       if(p4>p4max) p4max=p4;
       likeTot+=like;
       if(like>likeMax){
	 likeMax=like;
	 p1Best=p1;
	 p2Best=p2;
	 p3Best=p3;
	 p4Best=p4;
       }

    }
    fin.close();

    nentry=(int)(npoint/ndata);
    nz=nentry-nignore;

    cout<<npoint<<"\t"<<ndata<<"\t"<<npoint/ndata<<endl;
    cout<<p1min<<"\t"<<p1max<<"\t"<<p2min<<"\t"<<p2max<<"\t"<<p3min<<"\t"<<p3max<<"\t"<<p4min<<"\t"<<p4max<<endl;
       cout<<"Best p: "<<p0<<"\t"<<p1Best<<"\t"<<p2Best<<"\t"<<p3Best<<"\t"<<p4Best<<endl;

       pvol=p1max-p1min;
       if(fabs(p2max-p2min)>1.0e-4) pvol*=p2max-p2min;
       if(fabs(p3max-p3min)>1.0e-4) pvol*=p3max-p3min;
       if(fabs(p4max-p4min)>1.0e-4) pvol*=p4max-p4min;
       evidence=likeTot/(double)(ndata)*pvol/pvol;
       cout<<"Evidence: "<<evidence<<endl;
       cout<<"Max like: "<<likeMax<<endl;
       cout<<"Best prob: "<<likeMax/likeTot<<endl;
       cout<<"Best prob renorm: "<<likeMax/likeTot*(double)(ndata)<<endl;       
       ////////////////////////////
       //calculate marginalised pdfs using the above sampling points
       ////////////////////////////

    //calculate p1 pdf
    for(j=1;j<=nparam;j++){
      cout<<"p"<<j<<" pdf"<<endl;
    nbins=20;
    nbins=40;
    pdf.clear();
    for(i=0;i<nbins;i++) pdf.push_back(0.0);

    if(j==1){
      xmin=p1min;
      xmax=p1max;
      pv=p1v;
    }else if(j==2){
      xmin=p2min;
      xmax=p2max;
      pv=p2v;
    }else if(j==3){
      xmin=p3min;
      xmax=p3max;
      pv=p3v;
    }else if(j==4){
      xmin=p4min;
      xmax=p4max;
      pv=p4v;
    }

    dx=(xmax-xmin)/(double)(nbins-1);

    for(i=0;i<ndata;i++){
      ibin=floor((pv[i]-xmin)/dx);
       pdf[ibin]+=likev[i]/likeTot;
    }

    cout<<pdf.size()<<endl;

    wnorm=0.0;
    for(i=0;i<nbins;i++){
      wnorm+=pdf[i]*dx;
    }

    if(j==1) files=data_dir+"p1_pdf.dat";
    if(j==2) files=data_dir+"p2_pdf.dat";
    if(j==3) files=data_dir+"p3_pdf.dat";
    if(j==4) files=data_dir+"p4_pdf.dat";

    fout.open(files.c_str());
    for(i=0;i<nbins;i++){
       fout<<xmin+i*dx<<"\t"<<pdf[i]/wnorm<<endl;
    }
    fout.close();

    }
    ///////////////////////////////////
    // Calculate tau cmb distribution
      if(tau_flag==1){

	cout<<"tau"<<j<<" pdf"<<endl;
	nbins=20;
	nbins=40;

	xmin=0.0;
	xmax=0.5;    
	dx=(xmax-xmin)/(double)(nbins-1);
	iz=nignore;

	files=data_dir+file_like;
	xv.clear();
	fin.open(files.c_str()); 
	while(!fin.eof()){
	  for(i=1;i<nignore;i++) fin>>trash;
	  fin>>tau;
	  xv.push_back(tau);
	  getline(fin,str);
	}
	fin.close();

       pdf.clear();
       pdf.reserve(nbins);
       for(i=0;i<nbins;i++) pdf[i]=0.0;

	for(i=0;i<ndata;i++){
	  ibin=floor((xv[i]-xmin)/dx);
	  pdf[ibin]+=likev[i]/likeTot;
	}

	cout<<pdf.size()<<endl;

	wnorm=0.0;
	for(i=0;i<nbins;i++){
	  wnorm+=pdf[i]*dx;
	}

	files=data_dir+"tau_pdf.dat";

	fout.open(files.c_str());
	for(i=0;i<nbins;i++){
	  fout<<xmin+i*dx<<"\t"<<pdf[i]/wnorm<<endl;
	}
	fout.close();

      }


    //////////////////////////
    
   
       //calculate xi pdf at given redshift
    
    files=data_dir+"xi_contour.dat";
    fout_xi.open(files.c_str());


       z=7.0;
       z=6.25;
    while(z<18.0){
       cout<<z<<endl;


       iz=floor((z-zmin)/dz)+nignore;

       files=data_dir+file_like;
       xv.clear();
       fin.open(files.c_str()); 
       while(!fin.eof()){
          for(i=1;i<iz;i++) fin>>trash;
          fin>>xi;
          xv.push_back(xi);
          getline(fin,str);
       }
       fin.close();

       nbins=100;
       pdf.clear();
       pdf.reserve(nbins);
       for(i=0;i<nbins;i++) pdf[i]=0.0;
       xmin=0.0;
       xmax=1.0;
       dx=(xmax-xmin)/(nbins-1);

       for(i=0;i<ndata;i++){
          ibin=floor((xv[i]-xmin)/dx);
          pdf[ibin]+=likev[i]/likeTot/dx;
       }
    
       files=data_dir+"xi_pdf_z"+numToString(z)+".dat";
       fout.open(files.c_str());
       pdftot=0.0;
       for(i=0;i<nbins;i++){
	 pdftot+=pdf[i]*dx;
	 fout<<xmin+i*dx<<"\t"<<pdf[i]<<"\t"<<pdftot<<endl;

	 //cumulative pdf
	 //pTB=rz.getPowerFZH(z,k,xmin+i*dx);  //fix xi=fcoll*zeta
	 fout_xi<<z<<"\t"<<xmin+i*dx<<"\t"<<pdftot<<"\t"<<tocm.tBrightSat(z,xmin+i*dx)<<"\t"<<tocm.tBrightSat(z,xmin+i*dx)/tocm.tSky(z)<<endl;
       }
       fout.close();

       //z+=1.0;
       z+=dz;
    }
    fout_xi.close();
    
    //calculate redshift distribution of xi    
       xi=0.1;
       while(xi<0.99){
          cout<<xi<<endl;

          nbins=60;
          pdf.clear();
          pdf.reserve(nbins);
          for(i=0;i<nbins;i++) pdf[i]=0.0;
          xmin=6.0;
	  if(nocut_flag==1) xmin=4.0;
          xmax=20.0;
          dx=(xmax-xmin)/(nbins-1);

          files=data_dir+file_like;
          fin.open(files.c_str()); 
          for(i=0;i<ndata;i++){
	    z=0.0;
             xinew=1.1;
             for(j=1;j<=nignore;j++) fin>>trash;
             for(j=0;j<nz;j++){
                xiold=xinew;
                fin>>xinew;
                //linearly interpolate between bins
                if(xi<xiold && xi>xinew){
		  z=zmin+(j-1)*dz+dz*(xi-xiold)/(xinew-xiold);
		  //cout<<i<<"\t"<<z<<"\t"<<likev[i]<<endl;
		}
             }
             getline(fin,str);
             ibin=floor((z-xmin)/dx);
	     if(ibin>=0 && ibin<nbins) pdf[ibin]+=likev[i]/likeTot/dx;       
          }
          fin.close();

          wnorm=0.0;
          for(i=0;i<nbins;i++) wnorm+=pdf[i]*dx;
          //cout<<xi<<"\t"<<wnorm<<endl;

          files=data_dir+"z_pdf_xi"+numToString(xi)+".dat";
          fout.open(files.c_str());
          for(i=0;i<nbins;i++){
             fout<<xmin+i*dx<<"\t"<<pdf[i]<<endl;
          }
          fout.close();

          xi+=0.1;
       }
   

       //repeat for the final reionization bin p(xi>0.95)
       //note this is probability in each bin P(xi>0.95|z) 
       //not the pdf i.e. dP/dz as in the above
       xi=0.95;
          cout<<xi<<endl;

          nbins=50;
          pdf.clear();
          pdf.reserve(nbins);
          for(i=0;i<nbins;i++) pdf[i]=0.0;
          xmin=6.0;
	  if(nocut_flag==1) xmin=4.0;
          xmax=20.0;
	  dx=0.25;

          files=data_dir+file_like;
          fin.open(files.c_str()); 
          for(i=0;i<ndata;i++){
	    z=0.0;
             xinew=1.1;
             for(j=1;j<=nignore;j++) fin>>trash;
             for(j=0;j<nz;j++){
                xiold=xinew;
                fin>>xinew;
                //if xi > xi threshold add in
                if(xinew>xi){
		  z=zmin+j*dx;
		  //cout<<i<<"\t"<<z<<"\t"<<likev[i]<<endl;
		  ibin=floor((z-xmin)/dx);
		  if(ibin>=0 && ibin<nbins) pdf[ibin]+=likev[i]/likeTot; 
		}
             }
	     getline(fin,str);
          }
          fin.close();

          wnorm=0.0;
          for(i=0;i<nbins;i++) wnorm+=pdf[i]*dx;
          //cout<<xi<<"\t"<<wnorm<<endl;

          files=data_dir+"z_pdf_xi0.95.dat";
          fout.open(files.c_str());
          for(i=0;i<nbins;i++){
             fout<<xmin+i*dx<<"\t"<<pdf[i]<<endl;
          }
          fout.close();

       /////////////////////////////////////////////////
	 // Calculate key to map ionization history to fluctuations
	 ///////////////////////////////////////////////
	 /*
	 double pTB, k(0.1);
	 //pTB=rz.getPowerFZH(z,k,xmin+i*dx);  //fix xi=fcoll*zeta
	 files=data_dir+"fluc_key.dat";
	 fout.open(files.c_str());
	 
	 nbins=20;
	 dx=1.0/(double)(nbins-1);

	   fout<<z;
	   for(i=0;i<nbins;i++) fout<<"\t"<<i*dx;
	   fout<<endl;

	 z=18.0;
	 while(z>6.0){
	   xidd(0.0,z,&c,1);
	   fout<<z;
	   for(i=0;i<nbins;i++){
	     xi=i*dx;
	     pTB=rz.getPowerFZH(z,k,xi);  //fix xi=fcoll*zeta
	     fout<<"\t"<<pTB;
	     cout<<"\t"<<pTB<<endl;
	   }
	   fout<<endl;
	   z-=1.0;
	 }
	 fout.close();
	 */
       /////////////////////////////////////////////////
       //calculate bestfit history
       /////////////////////////////////////////////////
       double zstep;
       double Nion,gamma,zeta,clump,mfp;
       double Nion0(1e51);

       param.clear();
       param.push_back(p0);
       param.push_back(p1Best);
       param.push_back(p2Best);
       param.push_back(p3Best);
       param.push_back(p4Best);
       Ion.setZetaParam(0.0,param,1);

       files=data_dir+"best_history.dat";
       fout.open(files.c_str());

       cout<<"Best p: "<<p0<<"\t"<<p1Best<<"\t"<<p2Best<<"\t"<<p3Best<<"\t"<<p4Best<<endl;
       cout<<"tau="<<Ion.getTauCMB()<<endl;

       fout<<"#"<<p0<<"\t"<<p1Best<<"\t"<<p2Best<<"\t"<<p3Best<<"\t"<<p4Best<<endl;
       fout<<"# like_best like_tot best_prob evidence prior"<<likeMax<<"\t"<<likeTot<<"\t"<<likeMax/likeTot<<"\t"<<evidence<<"\t"<<1.0/pvol<<endl;
       fout<<"# tau="<<Ion.getTauCMB()<<endl;
       fout<<"# z Nion  xi  gamma  zeta  mfp  clump"<<endl;
       
       z=25.0;
       z=40.0;
       zstep=0.5;
       while(z>2.5){
	 Nion=Ion.getNion(z)*MPC*MPC*MPC;
	 xi=Ion.getXI(z);
	 gamma=lyf.getGammaFromNdot(z,Nion);
	 zeta=Ion.getZeta(z);
	 mfp=lyf.mfpFromGamma(z,gamma);
	 clump=Ion.getClumping(z,xi);

	 cout<<z<<"\t"<<Nion/Nion0<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<mfp<<"\t"<<clump<<endl;
	 fout<<z<<"\t"<<Nion/Nion0<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<mfp<<"\t"<<clump<<endl;	
	 z-=zstep;
       }

       fout.close();

///////////////////////////////////////////////////////////////////////
  return 0;
}




