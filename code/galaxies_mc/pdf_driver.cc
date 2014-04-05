
 /* driver.cc
 * program for mimicing the results of BL04
 * Call by 
 *   driver.x -z6 > (output filename)
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
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
#include "galaxies.h"

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

  double zin_arg(4.0);
  int tin(2);
  int i,j;
  char *file;
  ofstream fout;
  ifstream fin;
  int xin(1);
  int lyaxray_in(0);
  int popflag_in(0);
  int extern_in(0);

     int ndot_flag(0), cmb_flag(0);
     int param_flag(0);
     int lls_flag(1);
     int nthread(0);
     int lamLLS_flag(1);
     int case_flag(0);
     int m_flag(0);

  // Handle arguments 
  while ((argc>1) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'z':
      zin_arg=atof(&argv[1][2]);
      break;
    case 't':
      tin=atoi(&argv[1][2]);
      break;
    case 'e':
      extern_in=atoi(&argv[1][2]);
      break;
    case 'x':
      xin=atoi(&argv[1][2]);
      break;
    case 'l':
      lls_flag=atoi(&argv[1][2]);
      break;
    case 'm':
      lamLLS_flag=atoi(&argv[1][2]);
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
    case 'a':
      case_flag=atoi(&argv[1][2]); 
      break;
    case 'q':
      m_flag=atoi(&argv[1][2]); 
      break;
    default:
      cerr << "Bad option" <<argv[1] <<"\n";
    }
    --argc;
    ++argv;
  }

  
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
  Galaxies Gal;

  lyf.setLambdaLLS(lamLLS_flag);

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
  Gal.initGalaxies(&c,&a);



////////////////////////////////////////////////////////////////////////
// Convert gaussian dist x into new function y
///////////////////////////////////////////////////////////////////////
    /*
     string data_dir;
     string files;
     double z;
     double tau(0.09), sig_tau(0.017);
     long npoints;
     Ionization Ion(&c,&a);
     long idums;
     int *inxbin, nxbin, ixbin;
     double x, xmax, xmin, xbin;
     int *inybin, nybin, iybin;
     double y, ymax, ymin, ybin;
     

     npoints=100000;

     //initialise random number generator
     idums= -1;
     dran2(&idums);

     
     //x bins for underlying distribution
     xmin=tau-4*sig_tau;
     xmax=tau+4*sig_tau;
     nxbin=40;
     inxbin=ivector(0,nxbin-1);
     xbin=(xmax-xmin)/(double)(nxbin-1);

     //y bins for derived distribution
     ymin=0.0;
     ymax=8.0;
     nybin=40;
     inybin=ivector(0,nybin-1);
     ybin=(ymax-ymin)/(double)(nybin-1);

     //calculate distribution
     for(i=1;i<=npoints;i++){
       x=gaus(&idums,tau,sig_tau);
       ixbin=floor((x-xmin)/xbin);
       inxbin[ixbin]++;

       y=x/2.0;
       iybin=floor((y-ymin)/ybin);
       inybin[iybin]++;      
     }

     files="xpdf.dat";
     fout.open(files.c_str());
     for(i=1;i<=nxbin;i++) fout<<xmin+(i-1)*xbin+xbin/2.0<<"\t"<<inxbin[i-1]<<"\t"<<(double)(inxbin[i-1])/npoints/xbin<<endl;
     fout.close();
  
     files="ypdf.dat";
     fout.open(files.c_str());
     for(i=1;i<=nybin;i++) fout<<ymin+(i-1)*ybin+ybin/2.0<<"\t"<<inybin[i-1]<<"\t"<<(double)(inybin[i-1])/npoints/ybin<<endl;
     fout.close();
    */
////////////////////////////////////////////////////////////////////////
// Map gaussian error on tau_eff into Gamma and Ndot
///////////////////////////////////////////////////////////////////////
    /*
     string data_dir;
     string files;
     double zin(4.0), tau(0.805), sig_tau(0.067);  //bh values
     long npoints;
     Ionization Ion(&c,&a);
     long idums;
     int *inxbin, nxbin, ixbin;
     double x, xmax, xmin, xbin;
     int *inybin, nybin, iybin;
     double y, ymax, ymin, ybin;
     int *inzbin, nzbin, izbin;
     double z, zmax, zmin, zbin;

     int extern_flag(2);	
     int Tunif_flag(0);

     //external parameters
     double T0in(1.0e4), sig_T0(2.0e3), T0(T0in);
     double T0min(1.0e4), T0max(3.0e4);
     double beta_min(0.0), beta_max(0.6), beta(0.3);
     double betain, sig_beta;
     double aSmin(1.0), aSmax(3.0);
     double aS(3.0), aB(3.0);
     double aSin, sig_aS;
     double sig_nlls(0.2), nllsin(1.0), nlls(nllsin);
     double cindexin(1.5), sig_cindex(0.2), cindex(cindexin);

     string tag;
     stringstream ss;
     data_dir="";

     zin=6.0;
     zin=zin_arg;
     extern_flag=extern_in;
     npoints=50000;
     npoints=10000;

     if(case_flag==1){
       // model A
       data_dir="./pdf_data/modelA/";
       npoints=1000000;
       Tunif_flag=1;
       extern_flag=4;
     }else if(case_flag==2){
       //model B
       data_dir="./pdf_data/modelB/";
       npoints=1000000;
       Tunif_flag=1;
       extern_flag=2;
     }else if(case_flag==3){
       //model C - improve everything and see what happens
       data_dir="./pdf_data/modelC/";
       npoints=1000000;
       Tunif_flag=0;
       extern_flag=4;
     }

     //redshift dependent tau_eff data
     if(zin>4.1){
       tau=2.07;
       sig_tau=0.27;

       T0min=0.5e4;
       T0max=2.5e4;
     }
     if(zin>5.1){
       tau=5.50;
       sig_tau=10.00;  //really upperlimit here
       sig_tau=8.0;
       sig_tau=15.0;
       T0min=0.5e4;
       T0max=2.5e4;
     }

     if(case_flag==3){
       T0in=2.0e4;
       sig_T0=0.25e4;
       
       betain=0.3;
       sig_beta=0.15;

       aSin=1.5;
       sig_aS=0.25;

       sig_nlls=0.1;
       sig_tau/=2.0;
     }

     //file tag
     ss<<"_astro_z"<<zin;
     ss>>tag;

     //cout<<zin<<"\t"<<tag<<endl;

     //initialise random number generator
     idums= -1;
     dran2(&idums);
     
     //x bins for underlying distribution
     xmin=tau-4*sig_tau;
     xmax=tau+4*sig_tau;
     if(zin>5.1){ xmin=5.4; xmax=15.0;}
     nxbin=160;
     inxbin=ivector(0,nxbin-1);
     xbin=(xmax-xmin)/(double)(nxbin-1);

     //y bins for derived distribution
     ymin=0.2;
     ymax=1.9;
     if(zin>4.1){ ymin=0.0; ymax=1.9;}
     if(zin>5.1){ ymin=0.0; ymax=1.0; ymax=0.7;}
     nybin=160;
     inybin=ivector(0,nybin-1);
     ybin=(ymax-ymin)/(double)(nybin-1);

     //z bins for derived distribution
     zmin=0.0;
     zmax=2.2;
    if(zin>4.1){ zmin=0.1; zmax=2.8;}
    if(zin>5.1){ zmin=0.1; zmax=1.5;}
     nzbin=160;
     inzbin=ivector(0,nzbin-1);
     zbin=(zmax-zmin)/(double)(nzbin-1);


       files=data_dir+"log_pdf"+tag+".dat";
       fout.open(files.c_str());
       if(extern_flag>0) {if(Tunif_flag==0){fout<<"gaus T0\t"<<T0in<<"\t"<<sig_T0<<endl;
	 }else{ fout<<"unif T0\t"<<T0min<<"\t"<<T0max<<endl;}}
       if(extern_flag>1)fout<<"unif beta\t"<<beta_min<<"\t"<<beta_max<<endl;
       if(extern_flag>2)fout<<"unif aS\t"<<aSmin<<"\t"<<aSmax<<endl;
       if(extern_flag>3)fout<<"gaus nlls_norm\t"<<nlls<<"\t"<<sig_nlls<<endl;
       fout.close();

       //temp
       //ofstream fout2;


     //calculate distribution
     for(j=1;j<=npoints;j++){

       //first specify all parameters needed by code

       //external param
       if(extern_flag>0) T0=fmax(fmin(gaus(&idums,T0in,sig_T0),3.0e4),0.5e4);
       if(extern_flag>0 && Tunif_flag==1) T0=unif(&idums,T0min,T0max);
       if(extern_flag>1) beta=unif(&idums,beta_min,beta_max);
       if(extern_flag>2) {aS=unif(&idums,aSmin,aSmax); aB=aS;}
       if(extern_flag>3) nlls=gaus(&idums,nllsin, sig_nlls);
       if(extern_flag>4) cindex=gaus(&idums,cindexin, sig_cindex);

       if(extern_flag>1 && case_flag==3) beta=fmin(fmax(gaus(&idums,betain,sig_beta),beta_min),beta_max);
       if(extern_flag>1 && case_flag==3) {aS=fmax(fmin(gaus(&idums,aSin,sig_aS),aSmax),aSmin); aB=aS;}


       //cout<<aS<<endl;
       //nlls=gaus(&idums,3.3,0.6);

       //cosmology
       h=gaus(&idums,0.7,0.04);
       s8=gaus(&idums,0.8,0.05);
       om0=gaus(&idums,0.3,0.04);
       omb=gaus(&idums,0.046,0.0005);
       c.resetCosmology(om0,lam0,omb,h,s8,n,omNu);
       a.setCosmology(&c);
       lyf.initLyman(&c,&a);

       //lyf.setNormNLLS(1.0);
       lyf.setNormNLLS(nlls);
       //nlls=3.3;
       //lyf.setNormNLLS(nlls/lyf.dNLLSdz(2.4e4,y,zin));
       //cout<<y<<"\t"<<nlls<<"\t"<<lyf.dNLLSdz(2.4e4,y,zin)<<"\t"<<lyf.dNLLSdz(2.4e4,y,zin)/nlls<<endl;
       lyf.setTempDens(T0,beta);
       lyf.setColumnIndex(cindex);

       //now calculate observables
       //tau_eff
       x=gaus(&idums,tau,sig_tau);
       if(zin>5.1) x=unif(&idums,tau,sig_tau);

       ixbin=floor((x-xmin)/xbin);
       if(ixbin>nxbin-1) ixbin=nxbin-1;
       inxbin[ixbin]++;

       //fout2.open("temp_log.dat");
       //fout2<<x<<"\t"<<T0<<"\t"<<beta<<"\t"<<aS<<"\t"<<nlls<<"\t"<<cindex<<endl;     
       //fout2.close();

       //Gamma
       y=lyf.gammaFromTau(zin,x);
       iybin=floor((y-ymin)/ybin);
       if(iybin>nybin-1) iybin=nybin-1;
       inybin[iybin]++;   

       //fout2.open("temp_log.dat",ios::app);
       //fout2<<y<<endl;
       //fout2.close();

       //Ndot
       z=lyf.nionFromGamma(zin,y,aS,aB)/1.0e51;
       izbin=floor((z-zmin)/zbin);
       if(izbin>nzbin-1) izbin=nzbin-1;
       inzbin[izbin]++;  

       //fout2.open("temp_log.dat",ios::app);
       //fout2<<z<<endl;
       //fout2.close();

       if(j%1000==0){
	 cout<<zin<<"\t"<<j<<endl;
	 files=data_dir+"xpdf"+tag+".dat";
	 //cout<<tag<<"\t"<<files<<endl;
	 fout.open(files.c_str());
	 for(i=1;i<=nxbin;i++) fout<<xmin+(i-1)*xbin+xbin/2.0<<"\t"<<inxbin[i-1]<<"\t"<<(double)(inxbin[i-1])/j/xbin<<endl;
	 fout.close();
  
	 files=data_dir+"ypdf"+tag+".dat";
	 fout.open(files.c_str());
	 for(i=1;i<=nybin;i++) fout<<ymin+(i-1)*ybin+ybin/2.0<<"\t"<<inybin[i-1]<<"\t"<<(double)(inybin[i-1])/j/ybin<<endl;
	 fout.close();

	 files=data_dir+"zpdf"+tag+".dat";
	 fout.open(files.c_str());
	 for(i=1;i<=nzbin;i++) fout<<zmin+(i-1)*zbin+zbin/2.0<<"\t"<<inzbin[i-1]<<"\t"<<(double)(inzbin[i-1])/j/zbin<<endl;
	 fout.close();
       }

     
     }

     //finish up
     files=data_dir+"xpdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nxbin;i++) fout<<xmin+(i-1)*xbin+xbin/2.0<<"\t"<<inxbin[i-1]<<"\t"<<(double)(inxbin[i-1])/npoints/xbin<<endl;
     fout.close();
  
     files=data_dir+"ypdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nybin;i++) fout<<ymin+(i-1)*ybin+ybin/2.0<<"\t"<<inybin[i-1]<<"\t"<<(double)(inybin[i-1])/npoints/ybin<<endl;
     fout.close();

     files=data_dir+"zpdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nzbin;i++) fout<<zmin+(i-1)*zbin+zbin/2.0<<"\t"<<inzbin[i-1]<<"\t"<<(double)(inzbin[i-1])/npoints/zbin<<endl;
     fout.close();

     //summarise constraints by finding 1 and 2 sigma bands and maximum value
     
     double sigm2(0.0227501), sigp2(0.97725), sigm1(0.158655), sigp1(0.841345);
     double p2m, p2p, p1m, p1p, pmax;
     int imax;
     double pcum(0.0), pcumold(0.0);
     double *pbin;
     pbin=dvector(0,nzbin-1);

     files=data_dir+"errors"+tag+".dat";
     fout.open(files.c_str());

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nxbin;i++){
       pbin[i]=(double)(inxbin[i-1])/npoints/xbin;
       pcum+=pbin[i]*xbin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=xmin+(p2m-1)*xbin+xbin/2.0;
     p1m=xmin+(p1m-1)*xbin+xbin/2.0;
     p1p=xmin+(p1p-1)*xbin+xbin/2.0;
     p2p=xmin+(p2p-1)*xbin+xbin/2.0;
     pmax=xmin+(imax-1)*xbin+xbin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<xbin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;
     fout<<endl;

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nybin;i++){
       pbin[i]=(double)(inybin[i-1])/npoints/ybin;
       pcum+=pbin[i]*ybin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=ymin+(p2m-1)*ybin+ybin/2.0;
     p1m=ymin+(p1m-1)*ybin+ybin/2.0;
     p1p=ymin+(p1p-1)*ybin+ybin/2.0;
     p2p=ymin+(p2p-1)*ybin+ybin/2.0;
     pmax=ymin+(imax-1)*ybin+ybin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<ybin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;
     fout<<endl;

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nzbin;i++){
       pbin[i]=(double)(inzbin[i-1])/npoints/zbin;
       pcum+=pbin[i]*zbin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=zmin+(p2m-1)*zbin+zbin/2.0;
     p1m=zmin+(p1m-1)*zbin+zbin/2.0;
     p1p=zmin+(p1p-1)*zbin+zbin/2.0;
     p2p=zmin+(p2p-1)*zbin+zbin/2.0;
     pmax=zmin+(imax-1)*zbin+zbin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<zbin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;

     fout.close();

     free_ivector(inxbin,0,nxbin-1);
     free_ivector(inybin,0,nybin-1);
     free_ivector(inzbin,0,nzbin-1);

    */
////////////////////////////////////////////////////////////////////////
// Map uncertainties in galaxy LF into constraint on Ndot
///////////////////////////////////////////////////////////////////////

     string data_dir;
     string files;
     double zin(4.0), tau(0.805), sig_tau(0.067);  //bh values
     long npoints;
     Ionization Ion(&c,&a);
     long idums;
     int *inxbin, nxbin, ixbin;
     double x, xmax, xmin, xbin;
     int *inybin, nybin, iybin;
     double y, ymax, ymin, ybin;
     int *inzbin, nzbin, izbin;
     double z, zmax, zmin, zbin;

     double e25min(0.1), e25max(10.0);
     int extern_flag(0);	
     double Lmin, nuObs, lamObs(1500.0e-8);
     double mfp;
     double M0,M0in, sig_M0;

     //external parameters
     double phi0, sig_phi0, phi0in;
     double L0, sig_L0, L0in;
     double alpha, sig_alpha, alphain;
     double sig_fesc(1.0), fescin(0.02);
     double aS(3.0);
     double aSmin(1.0), aSmax(3.0);

     string tag;
     stringstream ss;
     data_dir="";

     zin=6.0;
     zin=zin_arg;
     extern_flag=extern_in;
     npoints=50000;
     //  npoints=10000;

     if(case_flag==1){
       // model A
       data_dir="./pdf_data/galA/";
       npoints=50000;
       extern_flag=0;
     }else if(case_flag==2){
       //model B
       data_dir="./pdf_data/galB/";
       npoints=50000;
       extern_flag=1;
     }else if(case_flag==3){
       //model C 
       data_dir="./pdf_data/galC/";
       npoints=50000;
       extern_flag=2;
     }

     //establish redshifts where I have info
     if(zin>3.5 && zin<4.5) z=3.8;
     if(zin>4.5 && zin<5.5) z=5.0;
     if(zin>5.5 && zin<6.5) z=5.9;
     if(zin>6.5 && zin<7.5) z=7.4;       
     if(zin>7.5) z=10.;

     //redshift dependent astrophysics parameters
      ss.clear();
       ss<<"./lf_z"<<z<<".dat";
       ss>>files;
       cout<<files<<endl;
       fin.open(files.c_str());
       fin>>zin>>lamObs;
       fin>>alphain>>sig_alpha;
       fin>>phi0in>>sig_phi0;
       fin>>M0in>>sig_M0;
       fin.close();
       z=zin;
       cout<<"reading in data for z="<<zin<<endl;
       cout<<alphain<<"\t"<<phi0in<<"\t"<<M0in<<endl;

     //file tag
       ss.clear();
     ss<<"_gal_z"<<zin;
     ss>>tag;

     //cout<<zin<<"\t"<<tag<<endl;

     //initialise random number generator
     idums= -1;
     dran2(&idums);
     
     //x bins for underlying distribution - epsilon25
     xmin=e25min;
     xmax=e25max;
     xmin=0.01;
     xmax=10.0;
     nxbin=160;
     inxbin=ivector(0,nxbin-1);
     xbin=(xmax-xmin)/(double)(nxbin-1);
     for(i=0;i<nxbin;i++) inxbin[i]=0;

     //y bins for derived distribution
     ymin=0.0;
     ymax=1.9;
     //if(zin>4.1){ ymin=0.0; ymax=1.9;}
     //if(zin>5.1){ ymin=0.0; ymax=1.0; ymax=0.7;}
     nybin=160;
     inybin=ivector(0,nybin-1);
     ybin=(ymax-ymin)/(double)(nybin-1);
     for(i=0;i<nybin;i++) inybin[i]=0;

     //z bins for derived distribution
     zmin=0.0;
     zmax=2.2;
     zmax=3.0;
     //if(zin>4.1){ zmin=0.1; zmax=2.8;}
     //if(zin>5.1){ zmin=0.1; zmax=1.5;}
     nzbin=160;
     inzbin=ivector(0,nzbin-1);
     zbin=(zmax-zmin)/(double)(nzbin-1);
     for(i=0;i<nzbin;i++) inzbin[i]=0;

     //output parameters used and ranges they can take
       files=data_dir+"log_pdf"+tag+".dat";
       fout.open(files.c_str());
       //put in data on parameters and values
       fout.close();

       //temp
       //ofstream fout2;

       //set lower limit of integrations
       //derived from LF data
       nuObs=SPEEDOFLIGHT_CGS/(lamObs);
       L0in=Gal.lumFromAbsMag(M0in,z)*nuObs;
       sig_L0=fabs(Gal.lumFromAbsMag(M0in+0.19,z)*nuObs-L0in);

       //setting lower limit is not obvious
       //but doesn't seem to make a large difference
       Lmin=Gal.lumFromAbsMag(-16.0,z)*nuObs;
       //Lmin=Gal.lumFromAbsMag(-19.0,z)*nuObs;
       Lmin=0.01*L0in;  //could set this for each run but not too big an effect

     //calculate distribution
     for(j=1;j<=npoints;j++){

       //first specify all parameters needed by code

       //cosmology
       h=gaus(&idums,0.7,0.04);
       s8=gaus(&idums,0.8,0.05);
       om0=gaus(&idums,0.3,0.04);
       omb=gaus(&idums,0.046,0.0005);
       c.resetCosmology(om0,lam0,omb,h,s8,n,omNu);
       a.setCosmology(&c);
       Gal.initGalaxies(&c,&a);

       //astrophysics parameters
       fesc=0.2;
       fesc=1.0;
       aS=3.0;
       phi0=phi0in;
       L0=L0in;
       alpha=alphain;
       mfp=40.0;

       phi0=fmax(gaus(&idums,phi0in,sig_phi0),0.0); //need to renormalise this
       //L0=gaus(&idums,L0in,sig_L0);
       M0=gaus(&idums,M0in,sig_M0);   // errors on M0 closer to Gaussian
       alpha=gaus(&idums,alphain,sig_alpha);
       L0=Gal.lumFromAbsMag(M0,zin)*nuObs;

       if(extern_flag>0) fesc=unif(&idums,fescin,sig_fesc);
       if(extern_flag>1) {aS=unif(&idums,aSmin,aSmax); fesc=0.2;}
       if(extern_flag>2) fesc=fmax(gaus(&idums,0.3,0.1),0.0);

       //cout<<Lmin<<"\t"<<M0<<"\t"<<L0<<"\t"<<phi0<<"\t"<<alpha<<"\t"<<fesc<<endl;
       //now calculate observables
       //e25
       x=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
       //cout<<x<<endl;

       ixbin=floor((x-xmin)/xbin);
       if(ixbin>nxbin-1) ixbin=nxbin-1;
       inxbin[ixbin]++;

       //Gamma
       y=Gal.gammaGal(x,aS,fesc,mfp,z);
       //cout<<y<<endl;

       iybin=floor((y-ymin)/ybin);
       if(iybin>nybin-1) iybin=nybin-1;
       inybin[iybin]++;   
       //cout<<y<<"\t"<<ymin<<"\t"<<ymax<<"\t"<<ybin<<"\t"<<iybin<<endl;

       //fout2.open("temp_log.dat",ios::app);
       //fout2<<y<<endl;
       //fout2.close();

       //Ndot
       z=Gal.ndotGal(x,aS,fesc)/1.0e51;
       //cout<<"z="<<z<<endl;
       izbin=floor((z-zmin)/zbin);
       if(izbin>nzbin-1) izbin=nzbin-1;
       inzbin[izbin]++;  
       //cout<<z<<"\t"<<zmin<<"\t"<<zmax<<"\t"<<zbin<<"\t"<<izbin<<endl;

       //fout2.open("temp_log.dat",ios::app);
       //fout2<<z<<endl;
       //fout2.close();

       if(j%1000==0){
	 cout<<zin<<"\t"<<j<<endl;
	 files=data_dir+"xpdf"+tag+".dat";
	 //cout<<tag<<"\t"<<files<<endl;
	 fout.open(files.c_str());
	 for(i=1;i<=nxbin;i++) fout<<xmin+(i-1)*xbin+xbin/2.0<<"\t"<<inxbin[i-1]<<"\t"<<(double)(inxbin[i-1])/j/xbin<<endl;
	 fout.close();
  
	 files=data_dir+"ypdf"+tag+".dat";
	 fout.open(files.c_str());
	 for(i=1;i<=nybin;i++) fout<<ymin+(i-1)*ybin+ybin/2.0<<"\t"<<inybin[i-1]<<"\t"<<(double)(inybin[i-1])/j/ybin<<endl;
	 fout.close();

	 files=data_dir+"zpdf"+tag+".dat";
	 fout.open(files.c_str());
	 for(i=1;i<=nzbin;i++) fout<<zmin+(i-1)*zbin+zbin/2.0<<"\t"<<inzbin[i-1]<<"\t"<<(double)(inzbin[i-1])/j/zbin<<endl;
	 fout.close();
       }

     
     }

     //finish up
     files=data_dir+"xpdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nxbin;i++) fout<<xmin+(i-1)*xbin+xbin/2.0<<"\t"<<inxbin[i-1]<<"\t"<<(double)(inxbin[i-1])/npoints/xbin<<endl;
     fout.close();
  
     files=data_dir+"ypdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nybin;i++) fout<<ymin+(i-1)*ybin+ybin/2.0<<"\t"<<inybin[i-1]<<"\t"<<(double)(inybin[i-1])/npoints/ybin<<endl;
     fout.close();

     files=data_dir+"zpdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nzbin;i++) fout<<zmin+(i-1)*zbin+zbin/2.0<<"\t"<<inzbin[i-1]<<"\t"<<(double)(inzbin[i-1])/npoints/zbin<<endl;
     fout.close();

     /////////////////////////////////////////////////////////////////
     //summarise constraints by finding 1 and 2 sigma bands and maximum value
     
     double sigm2(0.0227501), sigp2(0.97725), sigm1(0.158655), sigp1(0.841345);
     double p2m, p2p, p1m, p1p, pmax;
     int imax;
     double pcum(0.0), pcumold(0.0);
     double *pbin;
     pbin=dvector(0,nzbin-1);

     files=data_dir+"errors"+tag+".dat";
     fout.open(files.c_str());

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nxbin;i++){
       pbin[i]=(double)(inxbin[i-1])/npoints/xbin;
       pcum+=pbin[i]*xbin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=xmin+(p2m-1)*xbin+xbin/2.0;
     p1m=xmin+(p1m-1)*xbin+xbin/2.0;
     p1p=xmin+(p1p-1)*xbin+xbin/2.0;
     p2p=xmin+(p2p-1)*xbin+xbin/2.0;
     pmax=xmin+(imax-1)*xbin+xbin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<xbin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;
     fout<<endl;

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nybin;i++){
       pbin[i]=(double)(inybin[i-1])/npoints/ybin;
       pcum+=pbin[i]*ybin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=ymin+(p2m-1)*ybin+ybin/2.0;
     p1m=ymin+(p1m-1)*ybin+ybin/2.0;
     p1p=ymin+(p1p-1)*ybin+ybin/2.0;
     p2p=ymin+(p2p-1)*ybin+ybin/2.0;
     pmax=ymin+(imax-1)*ybin+ybin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<ybin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;
     fout<<endl;

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nzbin;i++){
       pbin[i]=(double)(inzbin[i-1])/npoints/zbin;
       pcum+=pbin[i]*zbin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=zmin+(p2m-1)*zbin+zbin/2.0;
     p1m=zmin+(p1m-1)*zbin+zbin/2.0;
     p1p=zmin+(p1p-1)*zbin+zbin/2.0;
     p2p=zmin+(p2p-1)*zbin+zbin/2.0;
     pmax=zmin+(imax-1)*zbin+zbin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<zbin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;

     fout.close();

     free_ivector(inxbin,0,nxbin-1);
     free_ivector(inybin,0,nybin-1);
     free_ivector(inzbin,0,nzbin-1);

////////////////////////////////////////////////////////////////////////
// Map uncertainties in QSO LF into constraint on Ndot
///////////////////////////////////////////////////////////////////////
/*
     string data_dir;
     string files;
     double zin(4.0), tau(0.805), sig_tau(0.067);  //bh values
     long npoints;
     Ionization Ion(&c,&a);
     long idums;
     int *inxbin, nxbin, ixbin;
     double x, xmax, xmin, xbin;
     int *inybin, nybin, iybin;
     double y, ymax, ymin, ybin;
     int *inzbin, nzbin, izbin;
     double z, zmax, zmin, zbin;

     double e25min(0.1), e25max(10.0);
     int extern_flag(0);	
     double Lmin, nuObs, lamObs(1450.0e-8);
     double mfp;
     double M0,M0in, sig_M0;

     //external parameters
     double phi0, sig_phi0, phi0in;
     double L0, sig_L0, L0in;
     double beta1, sig_beta1, beta1in;
     double beta2, sig_beta2, beta2in;
     double sig_fesc(1.0), fescin(0.02);
     double aS(3.0);
     double aSmin(1.0), aSmax(3.0);

     double gamma1, gamma1in, sig_gamma1;
     double gamma2, gamma2in, sig_gamma2;
     double phistar, phistarin, sig_phistar;

     string tag;
     stringstream ss;
     data_dir="";

     zin=6.0;
     zin=zin_arg;
     extern_flag=extern_in;
     npoints=50000;
       npoints=10000;

     if(case_flag==1){
       // model A
       data_dir="./pdf_data/qsoA/";
       npoints=50000;
       extern_flag=0;
     }else if(case_flag==2){
       //model B
       data_dir="./pdf_data/qsoB/";
       npoints=50000;
       extern_flag=1;
     }else if(case_flag==3){
       //model C 
       data_dir="./pdf_data/qsoC/";
       npoints=50000;
       extern_flag=2;
     }

     //establish redshifts where I have info
     if(zin>3.5 && zin<4.5) z=4.0;
     if(zin>4.5 && zin<5.5) z=5.0;
     if(zin>5.5 && zin<6.5) z=6.0;

     //redshift dependent astrophysics parameters
     //Bolton parameters
     if(m_flag==0){
      ss.clear();
       ss<<"./qso_lf_z"<<z<<".dat";
       ss>>files;
       cout<<files<<endl;
       fin.open(files.c_str());
       fin>>zin>>lamObs;
       fin>>phi0in>>sig_phi0;
       fin>>M0in>>sig_M0;
       fin>>beta1in>>sig_beta1;
       fin>>beta2in>>sig_beta2;
       fin.close();

     //Hopkins parameters
     }else{
      ss.clear();
       ss<<"./qso_hop_lf_z"<<z<<".dat";
       ss>>files;
       cout<<files<<endl;
       fin.open(files.c_str());
       fin>>zin>>lamObs;
       fin>>phistarin>>sig_phistar;
       fin>>L0in>>sig_L0;
       fin>>gamma1in>>sig_gamma1;
       fin>>gamma2in>>sig_gamma2;
       fin.close();

       phi0in=phistarin;
       sig_phi0=sig_phistar;

       beta1in= -(gamma1in+1.0);
       beta2in= -(gamma2in+1.0);
       sig_beta1=sig_gamma1;
       sig_beta2=sig_gamma2;
     }

       //check input was ok
       z=zin;
       cout<<"reading in data for z="<<zin<<endl;
       cout<<"\t"<<phi0in<<"\t"<<M0in<<"\t"<<beta1in<<"\t"<<beta2in<<endl;

     //file tag
       ss.clear();
     ss<<"_qso_z"<<zin;
     ss>>tag;

     //cout<<zin<<"\t"<<tag<<endl;

     //initialise random number generator
     idums= -1;
     dran2(&idums);
     
     //x bins for underlying distribution - epsilon25
     xmin=e25min;
     xmax=e25max;
     xmin=0.0001;
     xmax=12.0;
     nxbin=160;
     nxbin=500;
     inxbin=ivector(0,nxbin-1);
     xbin=(xmax-xmin)/(double)(nxbin-1);
     for(i=0;i<nxbin;i++) inxbin[i]=0;

     //y bins for derived distribution
     ymin=0.0;
     ymax=1.9;
     //if(zin>4.1){ ymin=0.0; ymax=1.9;}
     //if(zin>5.1){ ymin=0.0; ymax=1.0; ymax=0.7;}
     nybin=160;
     nybin=320;
     inybin=ivector(0,nybin-1);
     ybin=(ymax-ymin)/(double)(nybin-1);
     for(i=0;i<nybin;i++) inybin[i]=0;

     //z bins for derived distribution
     zmin=0.0;
     zmax=2.2;
     zmax=3.0;
     //if(zin>4.1){ zmin=0.1; zmax=2.8;}
     //if(zin>5.1){ zmin=0.1; zmax=1.5;}
     nzbin=160;
     nzbin=320;
     inzbin=ivector(0,nzbin-1);
     zbin=(zmax-zmin)/(double)(nzbin-1);
     for(i=0;i<nzbin;i++) inzbin[i]=0;

     //output parameters used and ranges they can take
       files=data_dir+"log_pdf"+tag+".dat";
       fout.open(files.c_str());
       //put in data on parameters and values
       fout.close();

       //temp
       //ofstream fout2;

       //set lower limit of integrations
       //derived from LF data
       nuObs=SPEEDOFLIGHT_CGS/(lamObs);

       //setting lower limit is not obvious
       //but doesn't seem to make a large difference
       Lmin=Gal.lumFromAbsMag(-22.0,z)*nuObs;
       //Lmin=0.01*L0in;

     //calculate distribution
     for(j=1;j<=npoints;j++){

       //first specify all parameters needed by code

       //cosmology
       h=gaus(&idums,0.7,0.04);
       s8=gaus(&idums,0.8,0.05);
       om0=gaus(&idums,0.3,0.04);
       omb=gaus(&idums,0.046,0.0005);
       c.resetCosmology(om0,lam0,omb,h,s8,n,omNu);
       a.setCosmology(&c);
       Gal.initGalaxies(&c,&a);

       //astrophysics parameters
       fesc=0.2;
       aS=3.0;
       phi0=phi0in;
       L0=L0in;
       M0=M0in;
       beta1=beta1in;
       beta2=beta2in;
       mfp=40.0;

       phi0=fmax(gaus(&idums,phi0in,sig_phi0),1.0e-9); //need to renormalise this
       M0=gaus(&idums,M0in,sig_M0);   // errors on M0 closer to Gaussian
       beta1=gaus(&idums,beta1in,sig_beta1);
       beta2=gaus(&idums,beta2in,sig_beta2);
       L0=Gal.lumFromAbsMag(M0,zin)*nuObs;
       if(m_flag==1){
          phi0=gaus(&idums,phi0in,sig_phi0);
          phi0=exp(phi0*log(10.0))/log(10.0);
          L0=gaus(&idums,L0in,sig_L0); 
          L0=exp(L0*log(10.0))*LSOLAR;
          Lmin=0.01*L0;
       }

       if(extern_flag>0) fesc=unif(&idums,fescin,sig_fesc);
       if(extern_flag>1) {aS=unif(&idums,aSmin,aSmax); fesc=0.2;}
       if(extern_flag>2) fesc=fmax(gaus(&idums,0.3,0.1),0.0);

       cout<<Lmin<<"\t"<<M0<<"\t"<<L0<<"\t"<<phi0<<"\t"<<beta1<<"\t"<<beta2<<"\t"<<fesc<<endl;
       //now calculate observables
       //e25
       x=Gal.emissivityQSO(Lmin,nuObs,z,phi0,L0,beta1,beta2);
       cout<<x<<endl;

       ixbin=floor((x-xmin)/xbin);
       if(ixbin>nxbin-1) ixbin=nxbin-1;
       inxbin[ixbin]++;

       //Gamma
       y=Gal.gammaQSO(x,aS,mfp,z);
       //cout<<y<<endl;

       iybin=floor((y-ymin)/ybin);
       if(iybin>nybin-1) iybin=nybin-1;
       inybin[iybin]++;   
       //cout<<y<<"\t"<<ymin<<"\t"<<ymax<<"\t"<<ybin<<"\t"<<iybin<<endl;

       //fout2.open("temp_log.dat",ios::app);
       //fout2<<y<<endl;
       //fout2.close();

       //Ndot
       z=Gal.ndotQSO(x,aS)/1.0e51;
       //cout<<"z="<<z<<endl;
       izbin=floor((z-zmin)/zbin);
       if(izbin>nzbin-1) izbin=nzbin-1;
       inzbin[izbin]++;  
       //cout<<z<<"\t"<<zmin<<"\t"<<zmax<<"\t"<<zbin<<"\t"<<izbin<<endl;

       //fout2.open("temp_log.dat",ios::app);
       //fout2<<z<<endl;
       //fout2.close();

       if(j%1000==0){
	 cout<<zin<<"\t"<<j<<endl;
	 files=data_dir+"xpdf"+tag+".dat";
	 //cout<<tag<<"\t"<<files<<endl;
	 fout.open(files.c_str());
	 for(i=1;i<=nxbin;i++) fout<<xmin+(i-1)*xbin+xbin/2.0<<"\t"<<inxbin[i-1]<<"\t"<<(double)(inxbin[i-1])/j/xbin<<endl;
	 fout.close();
  
	 files=data_dir+"ypdf"+tag+".dat";
	 fout.open(files.c_str());
	 for(i=1;i<=nybin;i++) fout<<ymin+(i-1)*ybin+ybin/2.0<<"\t"<<inybin[i-1]<<"\t"<<(double)(inybin[i-1])/j/ybin<<endl;
	 fout.close();

	 files=data_dir+"zpdf"+tag+".dat";
	 fout.open(files.c_str());
	 for(i=1;i<=nzbin;i++) fout<<zmin+(i-1)*zbin+zbin/2.0<<"\t"<<inzbin[i-1]<<"\t"<<(double)(inzbin[i-1])/j/zbin<<endl;
	 fout.close();
       }

     
     }

     //finish up
     files=data_dir+"xpdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nxbin;i++) fout<<xmin+(i-1)*xbin+xbin/2.0<<"\t"<<inxbin[i-1]<<"\t"<<(double)(inxbin[i-1])/npoints/xbin<<endl;
     fout.close();
  
     files=data_dir+"ypdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nybin;i++) fout<<ymin+(i-1)*ybin+ybin/2.0<<"\t"<<inybin[i-1]<<"\t"<<(double)(inybin[i-1])/npoints/ybin<<endl;
     fout.close();

     files=data_dir+"zpdf"+tag+".dat";
     fout.open(files.c_str());
     for(i=1;i<=nzbin;i++) fout<<zmin+(i-1)*zbin+zbin/2.0<<"\t"<<inzbin[i-1]<<"\t"<<(double)(inzbin[i-1])/npoints/zbin<<endl;
     fout.close();

     /////////////////////////////////////////////////////////////////
     //summarise constraints by finding 1 and 2 sigma bands and maximum value
     
     double sigm2(0.0227501), sigp2(0.97725), sigm1(0.158655), sigp1(0.841345);
     double p2m, p2p, p1m, p1p, pmax;
     int imax;
     double pcum(0.0), pcumold(0.0);
     double *pbin;
     pbin=dvector(0,nzbin-1);

     files=data_dir+"errors"+tag+".dat";
     fout.open(files.c_str());

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nxbin;i++){
       pbin[i]=(double)(inxbin[i-1])/npoints/xbin;
       pcum+=pbin[i]*xbin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=xmin+(p2m-1)*xbin+xbin/2.0;
     p1m=xmin+(p1m-1)*xbin+xbin/2.0;
     p1p=xmin+(p1p-1)*xbin+xbin/2.0;
     p2p=xmin+(p2p-1)*xbin+xbin/2.0;
     pmax=xmin+(imax-1)*xbin+xbin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<xbin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;
     fout<<endl;

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nybin;i++){
       pbin[i]=(double)(inybin[i-1])/npoints/ybin;
       pcum+=pbin[i]*ybin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=ymin+(p2m-1)*ybin+ybin/2.0;
     p1m=ymin+(p1m-1)*ybin+ybin/2.0;
     p1p=ymin+(p1p-1)*ybin+ybin/2.0;
     p2p=ymin+(p2p-1)*ybin+ybin/2.0;
     pmax=ymin+(imax-1)*ybin+ybin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<ybin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;
     fout<<endl;

     pcum=0.0;
     pcumold=0.0;
     pmax=0.0;
     for(i=1;i<=nzbin;i++){
       pbin[i]=(double)(inzbin[i-1])/npoints/zbin;
       pcum+=pbin[i]*zbin;
       if(pbin[i]>pmax){ pmax=pbin[i]; imax=i;}
       if(pcum>sigm2 && pcumold<sigm2) p2m=(double)(i-1)+(sigm2-pcumold)/(pcum-pcumold);
       if(pcum>sigm1 && pcumold<sigm1) p1m=(double)(i-1)+(sigm1-pcumold)/(pcum-pcumold);
       if(pcum>sigp1 && pcumold<sigp1) p1p=(double)(i-1)+(sigp1-pcumold)/(pcum-pcumold);
       if(pcum>sigp2 && pcumold<sigp2) p2p=(double)(i-1)+(sigp2-pcumold)/(pcum-pcumold);
       pcumold=pcum;
     }
     p2m=zmin+(p2m-1)*zbin+zbin/2.0;
     p1m=zmin+(p1m-1)*zbin+zbin/2.0;
     p1p=zmin+(p1p-1)*zbin+zbin/2.0;
     p2p=zmin+(p2p-1)*zbin+zbin/2.0;
     pmax=zmin+(imax-1)*zbin+zbin/2.0;
     fout<<p2m<<"\t"<<p1m<<"\t"<<pmax<<"\t"<<p1p<<"\t"<<p2p<<"\t"<<pcum<<"\t"<<zbin<<endl;
     fout<<pmax<<"\t"<<p1m-pmax<<"\t"<<p1p-pmax<<"\t"<<p2m-pmax<<"\t"<<p2p-pmax<<"\t"<<pcum<<endl;

     fout.close();

     free_ivector(inxbin,0,nxbin-1);
     free_ivector(inybin,0,nybin-1);
     free_ivector(inzbin,0,nzbin-1);

*/
///////////////////////////////////////////////////////////////////////
  return 0;
}




