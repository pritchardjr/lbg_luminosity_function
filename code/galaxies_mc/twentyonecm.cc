// twentyonecm.cc
// Contains function definitions needed for calculating 21cm class

#include <math.h>
#include <iostream>
#include <fstream>
#include "astrophysics.h"
#include "dcosmology.h"
#include "twentyonecm.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf.h>
#include "dnumrecipes.h"
#include <vector>

double xanorm(1.0);


using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".
//FOB (2006): Furlanetto, Oh, and Briggs (2006) Physics Reports
//BL2005 : Barkana and Loeb (2005) ApJ 626, 1

/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

TwentyOneCM::TwentyOneCM(Cosmology *c1, Astrophysics *a1)
{
  // cout <<"21cm constructor has been called." <<endl;
  c=c1;
  a=a1;
  //cout << "21cm Constructor called successfully." <<endl; 
  
}

TwentyOneCM::~TwentyOneCM()
{
  // cout <<"Atomic destructor has been called." <<endl;

}
/////////////////////////////////////////////////////////////////////////
// Member functions
/////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////
// Temperatures
/////////////////////////////////////////////////////////////////////////
double TwentyOneCM::tBright(double z)
{
  double tau,ts,tcmb,tbright;
  double *result;
  double tk,xi,xe,lya;

  result=dvector(1,3);
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];   
  xe=result[2];
  free_dvector(result,1,3);

  lya=a->lyaFlux(z);

  ts=tSpinGen(z,tk,xe,lya);
  tau=tau21CM(z,ts,xe);
  tcmb=TCMB*(1.0+z);
  tbright=tau*(ts-tcmb)/(1.0+z);

  return tbright;
}

double TwentyOneCM::tBrightGen(double z, double tk, double xi, double lya)
{
  double tau,ts,tcmb,tbright;

  ts=tSpinGen(z,tk,xi,lya);
  tau=tau21CM(z,ts,xi);
  tcmb=TCMB*(1.0+z);
  tbright=tau*(ts-tcmb)/(1.0+z);

  return tbright;
}

double TwentyOneCM::tSpin(double z)
{
  double tS,tG,tK;
  double xa,xc,ya,yc;
  double *result;
  double Xi,Xe,lya;

  result=dvector(1,3);
  a->getTIGM(z,result);
  tK=result[1];
  Xi=result[2];
  Xe=result[3];
  free_dvector(result,1,3);

  lya=a->lyaFlux(z);

  xa=getXAlpha(z,tK,Xe,lya);
  xc=getXColl(z,tK,Xe);

  tG=TCMB*(1.0+z);

  yc=xc*tG/tK;
  ya=xa*tG/tK;

  tS=tG+yc*tK+ya*tK;
  tS/=1+yc+ya;

  return tS;

}

double TwentyOneCM::tSpinGen(double z, double tK, double Xi, double lya)
{
  double tS,tG;
  double xa,xc,ya,yc;

  xa=getXAlpha(z,tK,Xi,lya);
  xc=getXColl(z,tK,Xi);

  tG=TCMB*(1.0+z);

  yc=xc*tG/tK;
  ya=xa*tG/tK;

  tS=tG+yc*tK+ya*tK;
  tS/=1.0+yc+ya;

  return tS;

}

//calculate the IGM temperature at redshift z
double TwentyOneCM::tKinetic(double z)
{
  double tk;
  double *result;
  result=dvector(1,3);
  a->getTIGM(z,result);
  tk=result[1];
  free_dvector(result,1,3);
  return tk;
}

//Calculate the ionization fraction at redshift z
double TwentyOneCM::getXFree(double z)
{
  double xe;
  double *result;
  result=dvector(1,3);
  a->getTIGM(z,result);
  xe=result[2];
  free_dvector(result,1,3);
  return xe;
}

//calculate the optical depth in the 21cm line
//Units 
double TwentyOneCM::tau21CM(double z, double tspin, double xi)
{
  double tau;

  tau=3.0*SPEEDOFLIGHT_CGS*PLANCKCONSTANT*A21CM*lambda21cm*lambda21cm;
  tau/=32.0*PI*BOLTZK*tspin*c->hubbleZ(z)*UNH*c->getH();
  tau*=(1.0-xi)*c->nh(z);

  return tau;
}

double TwentyOneCM::tBrightSat(double z, double xi)
{
  double tau,tcmb,tbright;
  double tk(1.0e4);

  tau=tau21CM(z,tk,xi);
  tcmb=TCMB*(1.0+z);
  tbright=tau*tk/(1.0+z);

  return tbright;
}

double TwentyOneCM::tSky(double z)
{
  double nu0(1.8e8);
  double T0(180.0);
  double Tsky,nu;

  nu=nu21cm/(1.0+z);
  Tsky=T0*pow(nu/nu0,-2.6);
  return Tsky;
}

///////////////////////////////////////////////////////////////////////
//coupling parameters
///////////////////////////////////////////////////////////////////////
//function to calculate the collisional coupling x_c
//slight problem here with definition of xi when He present
double TwentyOneCM::getXColl(double z, double tk, double xi)
{
  double xcoll(T21cm/A21CM);
  double rateEH, rateHH;
  xcoll /= TCMB*(1.0+z);

  rateEH=xi*kappaEH(tk);
  rateHH=(1.0-xi)*kappaHH(tk);

  xcoll *= c->nh(z)*(rateEH+rateHH);

  return xcoll;
}

double TwentyOneCM::getXCollEH(double z, double tk, double xi)
{
  double xcoll(T21cm/A21CM);
  double rateEH;
  xcoll /= TCMB*(1.0+z);

  rateEH=xi*kappaEH(tk);

  xcoll *= c->nh(z)*rateEH;

  return xcoll;
}

double TwentyOneCM::getXCollHH(double z, double tk, double xi)
{
  double xcoll(T21cm/A21CM);
  double rateHH;
  xcoll /= TCMB*(1.0+z);

  rateHH=(1.0-xi)*kappaHH(tk);

  xcoll *= c->nh(z)*rateHH;

  return xcoll;
}

// calculate the lyman alpha coupling, xalpha given a 
// proper Lya number flux lyaflux and a gas kinetic temperature tk.
//Units: tk        Kelvin
//       lyaflux   cm^-2 s^-1 Hz^-1 ster^-1
double TwentyOneCM::getXAlpha(double z, double tk, double xi, double lyaflux)
{
  double xalpha(16.0*PI*PI/27.0);
  double tcmb;
  
  tcmb=TCMB*(1.0+z);

  xalpha*=T21cm*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS*FALPHA;
  xalpha/=A21CM*tcmb*ELECTRONMASS*SPEEDOFLIGHT_CGS;

  xalpha*=getSAlpha(z,tk,xi);  //ME-Chen correction factor
  xalpha*=lyaflux;

  return xalpha;
}

//calculate lyman alpha correction factor
//I'll use fit from Chuzhoy and Shapiro
//
// NEED TO MAKE THIS MORE RIGOROUS - ie follow hirata
double TwentyOneCM::getSAlpha(double z, double tk, double xi)
{
  double salpha,alpha,tau;
  //  double a,gamma,eta;
  //double nualpha;

  //  deltanu=sqrt(2.0*BOLTZK*tk/ELECTRONMASS/SPEEDOFLIGHT_CGS/SPEEDOFLIGHT_CGS)*nualpha;
  //a=life/(4.0*PI*deltanu);

  //alpha=eta*pow(3.0*a/2.0/PI/gamma,0.33333);
  tau=3.0e5*(1.0-xi)*pow((1.0+z)/7.0,1.5);
  alpha=0.717*pow(tk,-2.0/3.0)*pow(tau/1.0e6,1.0/3.0);
  
  salpha=exp(-1.12*alpha);

  return salpha;
}

//////////////////////////////////////////////////////////////////////
// Calculate fluctuation coefficients
/////////////////////////////////////////////////////////////////////

// Calculate the coefficients in the expansion of \delta_{T_b}
// 1 : density
// 2 : neutral fraction
// 3 : lyman alpha coupling
// 4 : temperature 
// 5 : velocity gradient
//
// Reference: FOB (2006) Section 4
//
void TwentyOneCM::getBeta(double z, double tk, double xi, double lyaflux, double beta[])
{
  double xa,xc,xtot,xtott;
  double xcHH,xcEH;
  double Tcmb;

  Tcmb=TCMB*(1.0+z);
  xa=getXAlpha(z,tk,xi,lyaflux);
 
  xcEH= getXCollEH(z,tk,xi);
  xcHH= getXCollHH(z,tk,xi);
  xc=xcEH+xcHH;

  xtot=xa+xc;
  xtott=xtot*(1.0+xtot);

  beta[1]=1.0+xc/xtott;
  beta[2]=1.0+(xcHH-xcEH)/xtott;
  beta[3]=xa/xtott;
  beta[4]=Tcmb/(tk-Tcmb)+(xcEH*getDLogKappaEH(tk)+xcHH*getDLogKappaHH(tk))/xtott;
  beta[5]=1.0;

}

// Calculate the coefficients in the expansion of \delta T_b
// Note the difference from getBeta.  This essentially calculates 
// beta T_b to avoid divide by zero errors when T_b=0
// 1 : density
// 2 : neutral fraction
// 3 : lyman alpha coupling
// 4 : temperature 
// 5 : velocity gradient
//
// Reference: FOB (2006) Section 4
//
void TwentyOneCM::getBetaTb(double z, double tk, double xi, double lyaflux, double betaT[])
{
  double xa,xc,xtot,xtott;
  double xcHH,xcEH;
  double Tcmb,Tb,Ts,tau;
  double Tbmod;

  Tcmb=TCMB*(1.0+z);
  xa=getXAlpha(z,tk,xi,lyaflux);

 
  xcEH= getXCollEH(z,tk,xi);
  xcHH= getXCollHH(z,tk,xi);
  xc=xcEH+xcHH;

  xtot=xa+xc;
  xtott=xtot*(1.0+xtot);

  Tb=tBrightGen(z,tk,xi,lyaflux);

  betaT[1]=(1.0+xc/xtott)*Tb;
  betaT[2]=(1.0+(xcHH-xcEH)/xtott)*Tb;
  betaT[3]=(xa/xtott)*Tb;

  Ts=tSpinGen(z,tk,xi,lyaflux);
  tau=tau21CM(z,Ts,xi);
  Tbmod=tau/(1.0+z)*Ts*xtot/(1.0+xtot);  //Tb lacking (1-tg/tk) part

  betaT[4]=tau*Ts*xtot/(1.0+xtot)*Tcmb/tk/(1.0+z);//This term not zero even if Tb=0
  betaT[4]+=Tb*(xcEH*getDLogKappaEH(tk)+xcHH*getDLogKappaHH(tk))/xtott;
  betaT[5]=Tb;

}

//calculate the adiabatic index of the gas
double TwentyOneCM::getGammaA(double z, double tk)
{
  double gammaA;

  gammaA=2.0/3.0;   //adiabatic case 

  return gammaA;
}


////////////////////////////////////////////////////////////////////////
//  Calculate collision rate using Zygelman data
///////////////////////////////////////////////////////////////////////
//C code which needs to be updated sometime: taken from LAST

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

static gsl_spline *ZYGspline;
static gsl_interp_accel * ZYGacc;

// read in data and prepare spline 
int TwentyOneCM::Zygelman_init(void)
{
  int i, n;
  FILE *fp;
  double *T, *K;
  fp=fopen("Zygelman.dat", "r");
  if(fp==NULL){
    printf("can not open Zygelman.dat\n");
    exit(1);
  }
  n=0;
  while (fscanf(fp, "%*lg\t  %*le")!=EOF)
    n++;
  rewind(fp);

  T=dvector(0, n-1);
  K=dvector(0, n-1);
  for(i=0; i<n; i++)
    fscanf(fp, "%lg\t  %le", &(T[i]),  &(K[i]));
  fclose(fp);

  ZYGspline=gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(ZYGspline, T, K, n);
  return n;
}

//Calculate H-H collision rate: kappa_{1-0)^{HH}
// for any temperature T<1000, give collision induced spin change rate
//units:  T  Kelvin
//        R  cm^3 s^-1
double TwentyOneCM::kappaHH(double T)
{
  double R;
  static int N_ZYG=0;

  if(N_ZYG==0){
    N_ZYG=Zygelman_init();
    ZYGacc=gsl_interp_accel_alloc();
  }
  if(T<0){
    fprintf(stderr, "Zygelman: T=%g <0\n", T);
    exit(1);
  }

  if(T<1.0){
    cout<<"V. low temp in kappaHH"<<endl;
    return R=gsl_spline_eval(ZYGspline, 1.0, ZYGacc);
  }
  if(T<1000)
    R=gsl_spline_eval(ZYGspline, T, ZYGacc);
  else
    //Artifically apply fit to temperatures above T=1000K
    R=8.652e-11*pow(T,0.2067)*exp(-87.56/T);
  return R;
}


//Calculate e-H collision rate: kappa_{1-0)^{eH}
//Taken from Kuhlen 0510814, but originally from Liszt (2001)
//Units kappa cm^3 s^-1
//       tk   Kelvin
double TwentyOneCM::kappaEH(double tk)
{
  double kappa,temp;

  if(tk<1.0){
    cout<<"V. low temp in kappaEH"<<endl;
    tk=1.0;
  }

  if(tk>1.0e4){
    //    cout<<"V high temp in kappaEH"<<endl;
    tk=1.0e4;
  }

  temp=-9.607;
  temp+=0.5*log10(tk)*exp(-pow(log10(tk),4.5)/1800.0);
  kappa=exp(temp*log(10.0));

  return kappa;
}

// Logarithmic derivative of kappaHH
double TwentyOneCM::getDLogKappaHH(double tkin)
{
  double gDLK;
  double step(tkin*0.05);
  double lKp,lKm,dK;
  step/=2.0;

  lKp=log(kappaHH(tkin+step));
  lKm=log(kappaHH(tkin-step));
  dK=(lKp-lKm)/(2.0*step);
  gDLK=dK*tkin;

  return gDLK;
}

// Logarithmic derivative of kappaEH
double TwentyOneCM::getDLogKappaEH(double tkin)
{
  double gDLK;
  double step(tkin*0.05);
  double lKp,lKm,dK;
  step/=2.0;

  lKp=log(kappaEH(tkin+step));
  lKm=log(kappaEH(tkin-step));
  dK=(lKp-lKm)/(2.0*step);
  gDLK=dK*tkin;

  return gDLK;
}



////////////////////////////////////////////////////////////////////
// 21CM Cutoff scales
////////////////////////////////////////////////////////////////////
// DONT REALLY UNDERSTAND THIS SCALE
// TAKEN FROM BL2005

double TwentyOneCM::getCutoffThermal(double z, double tk)
{
  return 8.3*sqrt(tk/100.0)/sqrt((1.0+z)/7)*1.0e-3;
}

// Pressure smoothing scale of the IGM
// Shapiro, Giroux, Babul (1994)
// Note rF=lambda_F/(2 PI)
// Units: rF   Mpc
double TwentyOneCM::getCutoffFilter(double z)
{
  double rF;
  double mF,rho;

  mF=a->filterMassFull(z);
  rho=4.0*PI/3.0*CRITDENMSOLMPC*c->getOm0hh();
  rF=pow(mF/rho,1.0/3.0)/PI;  //mass defined by (lambda_F/2)^3

  return rF;
}

////////////////////////////////////////////////////////////////////
//  FLUCTUATION POWER SPECTRUM FROM THE FIRST STARS
////////////////////////////////////////////////////////////////////
// Code to calculate the brightness temperature fluctuations 
// from the first stars
// Barkana and Loeb 2004

int nmarkt,nmarkr;

//calculate 21cm brightness fluctuations from lya variation
double TwentyOneCM::powerSpectrumLya(double z, double k, double store[])
{
  double *betaT, *result;
  double tk,xi,xe,lyaflux;
  double wk;
  double pk,pkfull,ps,psfull;
  double DZ,rFcut,rTcut,rF,rT;
  double tb,gammaA;

  result=dvector(1,3);
  betaT=dvector(1,5);

  //calculate thermodynamic quantities
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lyaflux=a->lyaFlux(z);
  tb=(1.0-xi)*tBrightGen(z,tk,xe,lyaflux);
  getBetaTb(z,tk,xi,lyaflux,betaT);
  gammaA=getGammaA(z,tk);

  //get window function
  wk=getWindowLya(z,k);

  //cutoffs and power spectra
  DZ=c->growthFac(z);
  rF=getCutoffFilter(z);
  rT=getCutoffThermal(z,tk);
  rFcut=pow(1.0+k*k*rF*rF,-2.0);
  rTcut=exp(-k*k*rT*rT);

  ps=c->powerSpectrum(k)*DZ*DZ*rFcut*rTcut;
  pk=2.0*ps;
  pk*=betaT[5]*(betaT[1]+betaT[4]*gammaA+wk*betaT[3]);

  psfull=k*k*k*2.0*betaT[5]*(betaT[1]+betaT[4]*gammaA)*ps/2.0/PI/PI; 
  pkfull=k*k*k*pk/2.0/PI/PI;

  //output full results for curiosity
  store[1]=wk;
  store[2]=2.0*betaT[5]*(betaT[1]+betaT[4]*gammaA)*ps;
  store[3]=pk;
  store[4]=psfull;
  store[5]=pkfull;

  free_dvector(result,1,3);
  free_dvector(betaT,1,5);

  return pkfull;
}


//Calculate the density contribution to the lya fluctuations
double TwentyOneCM::getWindowLya(double z, double k)
{
  double wk(0.0);
  double zmin;
  double zmax,zmax2;
  double n;

  double npoint(2.0);
  double rzmax,zuse;
  double period_max(100.0);
  double b,cc;
  double zcut(30.0);
 
  double z2;
  double x_alpha,x_coll,x_tot,x_tot_tilde;
  double tk, lyaflux,xi, *result;
  double edge(0.0);
  double tol(1.0e-4);
  double xacut(40.0);  //don't calculate once xa>xacut
  xacut=300.0;   //cross-terms mean need to go to higher xacut

  //  cout<<z<<"\t"<<k<<endl;
  if(z>ZSFR) return 0.0;  //no lya flux -> no fluctuations
                          //might be better to do this using x_alpha
  if(k>KCUT) return 0.0;
  if(z<a->getZReion()) return 0.0;  // fluctuations not relevant after reionization

  result=dvector(1,3);

  a->getTIGM(z,result);
  tk=result[1];
  xi=result[3];
  lyaflux=a->lyaFlux(z);
 
  free_dvector(result,1,3);

  x_alpha = getXAlpha(z,tk,xi,lyaflux);
  if(x_alpha>xacut) return 0.0; //strongly coupled fluctuations unimportant

  x_coll=getXColl(z,tk,xi);
  x_tot=x_alpha+x_coll;
  x_tot_tilde=x_tot*(1.0+x_tot);

  //Calculate w(k)
  wk=0.0;
  n=2.0;
  setWindowLyaK(z,z,k,tk,xi,c,a,this,1);
  setDJalphaDz(z,z,c,a,1);
  //Calculate redshift limit for this line
  zmax2=a->lymanZMax(z,n);
  
  // Sanity check limits of integration: calculate step size
  // in comoving distance, but implement in redshift.
  rzmax=period_max*2.0*PI/k;
  if(rzmax<c->relativeConformalR(z,zmax2)) zcut=c->zOfR(rzmax,z);
  //  if(rzmax>c->relativeConformalR(z,ZSFR)) zcut=ZSFR;
  if(fabs(zcut-z)<1.0e-6) zcut=2.0*z;
  
  // integrate using trapezium rule
  zmin=z+a->getZHII(z);
  n=lynmax;
  zmax=a->lymanZMax(z,n);

  zuse=zmin;
  z2=c->zOfR(2.0*PI/k/npoint,zuse);
  b=0;
  cc=0;
  //edge=1.0e-4;
  while(1){
    if(z2>zmax){
      if(fabs(n-2.0)<1.0e-4) b=1;  //finish on Lya
      //	z2=0.99999*zmax;
      z2=zmax-edge;
      n-=1.0;
      zmax=a->lymanZMax(z,n);
      cc=0; //switch to skip discontinuities
    }
    if(z2>zcut){
      //	cout<<"hi"<<endl;
      z2=zcut;
      b=1;
    }
    wk+=qromb(intWindowLyaK,zuse,z2,tol);
    nmarkr++;
    if(b==1) break;
    zuse=z2;
    //drdz=SPEEDOFLIGHT_CGS/c->hubbleZ(zuse)/H0_CMSMPC/c->getH();
    //z2+=(2.0*PI/k/npoint)/drdz;
    z2=c->zOfR(2.0*PI/k/npoint,zuse);
    if(cc==1){
      zuse+=2.0*edge;
      cc=0;
    }
    }
  
  wk/=x_alpha;
  cout << z<<"\t"<<k <<"\t"<<nmarkr<<"\t"<<nmarkt<< "\t" <<wk<<"\t"<<zcut<<"\t"<<x_alpha<<endl;
  nmarkt=0;
  nmarkr=0;
  
  return wk;
}

//function for integrating window
double setWindowLyaK(double z1, double zuse, double k1,double tk1,double xi1,Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag)
{
  static Cosmology *c;
  static Astrophysics *a;
  static TwentyOneCM *tocm;
  static double z;
  static double xi;
  static double k;
  static double tk;
  double dJdz;
  double wk,rz,x,oscfun,temp;
  
  if(iflag==1){
    c=c1;
    a=a1;
    tocm=tocm1;
    z=z1;
    xi=xi1;
    k=k1;
    tk=tk1; //gas kinetic temperature
    return 0.0;
  }

  if(fabs(zuse-z)<1.0e-7){
    return 0.0;
  }
  rz=SPEEDOFLIGHT_CGS*(c->confTime(z)-c->confTime(zuse))/H0_CMSMPC/c->getH();
  x=rz*k;
  oscfun = (1.0+tocm->splineBias(zuse))*sin(x)/x;
  oscfun -= 2.0*((3.0/pow(x,3.0)-1.0/x)*sin(x)-3.0*cos(x)/x/x)/3.0;
  oscfun*= c->growthFac(zuse)/c->growthFac(z);

  //get total lya differential flux
  // dJdz=tocm->getDJalphaDz(z,zuse);
  dJdz=getDJalphaDz(zuse);
  temp=tocm->getXAlpha(z,tk,xi,dJdz);
  wk=temp*oscfun;
 
  nmarkt++;

  return wk;
  
}

//
double intWindowLyaK(double zuse)
{
  return setWindowLyaK(1.0,zuse,1.0,1.0,0.0,NULL,NULL,NULL,0);
}
//////////////////////////////////////////////////////////////////////
// Poisson fluctuation code - Lya
//////////////////////////////////////////////////////////////////////

//Calculate the Poisson contribution to the lya fluctuations
//because I have to calculate the correlation function and then 
//fourier transform it, its not worth getting any single k value
//I just do everything at once and store the correlation function
//and the power spectrum in files.

double TwentyOneCM::poissonPowerLya(double z)
{
  double pC;
  int  nr(20),nk,nuse(20);
  double l,lmin,lmax,dl;
  double k,kmin,kmax,dk;
  double *lvec,*pcvec;
  double *kvec,*ftvec;
  double tk,xi,xe,lyaflux,*betaT,*result;
  double rT,ps,pscut,psfull,rTcut;
  double smallfac(32.0);
  int i;
  char *file;
  ofstream fout;

  cout<<"r_alpha="<<c->relativeConformalR(z,a->lymanZMax(z,2.0))<<endl;

  //Thermodynamics
  result=dvector(1,3);
  betaT=dvector(1,5);
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lyaflux=a->lyaFlux(z);
  getBetaTb(z,tk,xe,lyaflux,betaT);

  //Set limits of correlation function
  lmax=2.0*c->relativeConformalR(z,a->lymanZMax(z,2.0));
  lmin=1.0e-4;
  dl=exp(log(lmax/lmin)/(double)(nuse));
  //find number of data points taken
  i=0;
  l=lmin;
  while(l<lmax){
    i++;
    if(l>lmax/4.0){ //when near end tighten up steps
      l*=1.0+(dl-1.0)/smallfac;
    }else{
      l*=dl;
    }
  }
  nr=i;
  cout<<lmin<<"\t"<<lmax<<"\t"<<dl<<"\t"<<nr<<endl;
  lvec=dvector(1,nr);
  pcvec=dvector(1,nr);

  //Calculate correlation function
  cout<<"Calculating Lya correlation function"<<endl;
  file="poisson_pc_lya.dat";
  fout.open(file);
  i=0;
  l=lmin;
  while(l<=lmax){
    i++;
    pC=poissonCorrelationLya(l,z);
    lvec[i]=l;
    pcvec[i]=pC;

    fout <<l<<"\t"<<pC<<endl;
    cout<<z<<"\t"<<l<<"\t"<<pC<<endl;

    if(l>lmax/4.0){ //when near end tighten up steps
      l*=1.0+(dl-1.0)/smallfac;
    }else{
      l*=dl;
    }
  }
  fout.close();

  //Now perform Fourier Transform
  nk=2*nuse;
  kvec=dvector(1,nk);
  ftvec=dvector(1,nk);

  kmax=1.0e3;;
  kmin=1.0e-3;
  dk=exp(log(kmax/kmin)/(double)(nk-1));

  //set kvector
  k=kmin;
  for(i=1;i<=nk;i++){
    kvec[i]=k;
    k*=dk;
  }
  cout <<"doing FT"<<endl;
  FTTable(pcvec,lvec,nr,ftvec,kvec,nk);
  cout<<"FT done"<<endl;

  file="poisson_ft_lya.dat";
  fout.open(file);
  for(i=1;i<=nk;i++){
    fout<<kvec[i]<<"\t"<<ftvec[i]<<endl;
   }
  fout.close();

  cout<<"calculating Lya power spectrum"<<endl;
  //calculate the cutoff scales
  rT=8.3*sqrt(tk/100.0)/sqrt((1.0+z)/7)*1e-3;
  file="poisson_tb_lya.dat";
  fout.open(file);
  fout<<"#"<<rT<<endl;
  for(i=1;i<=nk;i++){
    k=kvec[i];
    ps=ftvec[i];
    rTcut=exp(-k*k*rT*rT);
    pscut=ps*rTcut;
    pscut*=betaT[3]*betaT[3]; 
    psfull=k*k*k*pscut/2.0/PI/PI;
    cout << k <<"\t"<<ps<<"\t"<<pscut<<"\t"<<psfull<<endl;
    fout << k <<"\t"<<ps<<"\t"<<pscut<<"\t"<<psfull<<endl;   
  }
  fout.close();
  cout<<"Lya power spectrum done"<<endl;

  free_dvector(result,1,3);
  free_dvector(betaT,1,5);
  free_dvector(lvec,1,nr);
  free_dvector(pcvec,1,nr);
  free_dvector(kvec,1,nk);
  free_dvector(ftvec,1,nk);

  return rT;
}



double TwentyOneCM::weightPLya(double z,double zp, double tk, double xi)
{
  double wP,dJdz;
  double dzdr;

  dzdr=c->hubbleZ(zp)*H0_CMSMPC*c->getH()/SPEEDOFLIGHT_CGS;
  dJdz=getDJalphaDz(zp);
  wP=getXAlpha(z,tk,xi,dJdz)*dzdr;
  wP/=4.0*PI*a->splineFColl(zp,c)* CRITDENMSOLMPC*c->getOm0hh();
  return wP;
}

double setPoissonCorrelationKernelLya(double rA, double theta, double z1,double l1, double tk1, double xi1, Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag)
{
  double pC,rB;
  double zpA,zpB;
  double rlimit(300.0);
  static Cosmology *c;
  static Astrophysics *a;
  static TwentyOneCM *tocm;
  static double z;
  static double l;
  static double tk;
  static double xi;
  if(iflag==1){
    c=c1;
    a=a1;
    z=z1;
    l=l1;
    tk=tk1;
    xi=xi1;
    tocm=tocm1;
    return 0.0;
  }

  //  rlimit=relativeConformalR(z,lymanZMax(z,2.0));

  //impose the boundary limit 
  if(fabs(rA)<1.0e-6) {
    cout <<"using min bound"<<endl;
   return 0.0; 
  }
  if(rA*cos(theta)> l/2.0+1.0e-8)  {
    cout << "rA*cos(theta)>l/2" <<endl;
    //return 0.0;
  }
 
  zpA=c->zOfR(rA,z);
  rB=sqrt(l*l+rA*rA-2.0*l*rA*cos(theta));
  zpB=c->zOfR(rB,z);

  //NEED TO MODIFT THIS IN SOME WAY
  // avoid problems if out of spline reach, at these distances one of
  // P(zA) or P(zB) will be zero as beyond max distance contributing to Lya.
  if((rA>rlimit)||(rB>rlimit)){
    //cout <<"applying lyman bound in setPoissonCorrelationKernal"<<endl;
    return 0.0;
  }

  pC=rA*rA*sin(theta);
  pC*=tocm->weightPLya(z,zpA,tk,xi)*tocm->weightPLya(z,zpB,tk,xi)/rA/rA/rB/rB;
  pC*=a->splineFColl(zpB,c)/a->splineFColl(zpA,c);
  pC*=tocm->splineMFColl(zpA);

  return pC;
}

//dummy function for 2D integration
double intPoissonLya(double rA, double theta)
{
  return setPoissonCorrelationKernelLya(rA,theta,0.0,0.0,0.0,0.0,NULL,NULL,NULL,0);
}

static double rAsave;
static double lsave;


double TwentyOneCM::poissonCorrelationLya(double l, double z)
{
  double pC(0.0),rmax;
  double tol(1.0e-3);
  double xalpha;
  double lrmin,lrmax;
  double ralpha;
  double tk,xi,xe,lyaflux,*result;
  lsave=l;

  //Basic thermodynamic quantities
  result=dvector(1,3);
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lyaflux=a->lyaFlux(z);
  xalpha=getXAlpha(z,tk,xe,lyaflux);

  setPoissonCorrelationKernelLya(1.0,1.0,z,l,tk,xe,c,a,this,1);
  setDJalphaDz(z,z,c,a,1);
  //limit integral to twice largest lyman alpha distance
  ralpha=c->relativeConformalR(z,a->lymanZMax(z,2.0));

  //NEED TO CONSIDER NEXT STEP CAREFULLY
  rmax=fmin(10.0*l,ralpha*0.99999);
  lrmin=log(1.0e-5);
  //limit integration to region where can get contribution
  if(l>ralpha)  lrmin=log(l-ralpha+1.0e-5);
  lrmax=log(rmax);
  
  //Calculate correlation function
  pC=qromb(intPoissonRLya,lrmin,log(rmax),tol);
  pC*=2.0*PI;
  pC*=2.0/xalpha/xalpha;

  cout <<lsave<<"\t"<<"\t"<<nmarkr<<"\t"<<nmarkt<<"\t"<<pC<<endl;
  nmarkr=0;
  nmarkt=0;

  free_dvector(result,1,3);

  return pC;
}

//integral over log r for Poisson fluctuations
double intPoissonRLya(double lrA)
{
  double iPR;
  double tol(1.0e-3);
  double thetamin(0.0);
  double thetamax(PI);
  double rA;
  double ralpha;
  
  rA=exp(lrA);
 //ralpha=2.0*max distance Lya can contribute
  ralpha=555.76/2.0;

  if(rA>ralpha) return 0.0;

  //restrict geometry
  if(rA>lsave/2.0){ //restrict to rA<rB
    thetamin=acos(lsave/2.0/rA);
  }
  if(rA+lsave>ralpha){ //restrict to overlap region
    thetamax=rA*rA+lsave*lsave-ralpha*ralpha;
    thetamax/=2.0*rA*lsave;
    if(fabs(thetamax)<=1.0) {
      thetamax=acos(thetamax);
    }else{
      return 0.0;
    }
    if(thetamax<=thetamin) return 0.0;
  }

  rAsave=rA;
  iPR=qromb(intPoissonThetaLya,thetamin,thetamax,tol);
  iPR*=rA;

  // cout <<lsave<<"\t"<<rA<<"\t"<<thetamark<<"\t"<<nmarkr<<"\t"<<nmarkt<<"\t"<<iPR<<endl;
  //nmarkt=0;
  nmarkr++;

  return iPR;
}

double intPoissonThetaLya(double theta)
{
  double iPT;
  nmarkt++;
  iPT= intPoissonLya(rAsave,theta);
  return iPT;
}

///////////////////////////////////////////////////////////////////////
// spline for the mass*collapse fraction 
//////////////////////////////////////////////////////////////////////
// Create a spline of the collapse function for an arbitrary mass function
// initialises spline on first run and runs quickly from then on
// xuse = redshift to evaluate fColl at.
//uses files to save overhead
double TwentyOneCM::splineMFColl(double xuse)
{
  static double *ly2;
  static double *lxs, *lys;
  static int nmaxs;
  static int iflag,icount;
  double yp1,ypn;
  double lxstep;
  int n;
  double yout,lxuse;
  int massfcn(MASSFCN); 
  char *file="./mfcoll_spline.dat";
  ifstream fin(file);
  ofstream fout;


  lxuse=log(xuse);

  //initialise spline on first pass
  if(iflag==0){
    if(icount>0){
      free_dvector(lxs,1,nmaxs);
      free_dvector(lys,1,nmaxs);
      free_dvector(ly2,1,nmaxs);
    }
    // store information in static variables
    nmaxs=150;
    //nmaxs=40;
    ly2=dvector(1,nmaxs);
    lxs=dvector(1,nmaxs);
    lys=dvector(1,nmaxs);
    cout <<"initialising mfcoll spline"<<endl;

    //if file doesn't exist then initialise by hand
    if(!fin){
      //////////
      //      double zmin(9.0);
      double zmin(ZLOW);
      double zmax(46.0);
      double lxmin,lxmax;
      lxmin=log(zmin);
      lxmax=log(zmax);
      lxstep=(lxmax-lxmin)/double(nmaxs-1);
      
      for(n=1;n<=nmaxs;n++){
	lxs[n]=lxmin+lxstep*(double)(n-1);
	lys[n]=log(fMColl(exp(lxs[n]),coolMass(c,exp(lxs[n])),c,massfcn));
      }
      // logarithmic derivatives at endpoints
      yp1=1.0e30;
      ypn=1.0e30;
      //initialise spline
      spline(lxs,lys,nmaxs,yp1,ypn,ly2);
      //store information for later reference
      fout.open(file);
      for(n=1;n<=nmaxs;n++){
	fout <<lxs[n]<<"\t"<<lys[n]<<"\t"<<ly2[n]<<endl;
      }
      fout.close();
    }else{
      for(n=1;n<=nmaxs;n++){
	fin >>lxs[n]>>lys[n]>>ly2[n];
      }
      fin.close();
    }
    splint(lxs,lys,ly2,nmaxs,lxuse,&yout);
    //icount++;
    iflag=1;
    return exp(yout);
  }

  // apply cubic spline
  if(iflag==1){
    if((lxuse<lxs[1])||(lxuse>lxs[nmaxs])) {
      cout <<"lx exceeds limits in splineMFColl"<<endl;
      return 0.0;
    }
    splint(lxs,lys,ly2,nmaxs,lxuse,&yout);
    return exp(yout);
  }

  if(iflag==2){
    free_dvector(lxs,1,nmaxs);
    free_dvector(lys,1,nmaxs);
    free_dvector(ly2,1,nmaxs);
    icount=0;
    return 0.0;
  }

  cout << "Error in splineStart"<<endl;
  return 0.0;
}

// Calculates fraction of mass in galaxies by directly integrating the
// mass function (using setFCollInt).  
//   z = redshift
//   mMin = minimum mass; if <0 use Tvir=10^4 K (optional)
//   massFcn = which mass function to use; 0=PS,1=ST,2=Jenkins (optional) 
//
double fMColl(double z, double mMin, Cosmology *c,int massFcn)
{
  double mMax;
  double ans,oldAns;

  setFMCollInt(0.0,c,z,massFcn,1);
  ans = 0.0;
  if (mMin < 0.0)
    mMin = coolMass(c,z);
  while (1) {
    mMax = mMin*10.0;
    oldAns = ans;
    ans += qromb(fMCollInt,mMin,mMax,1.0e-4);
    if (fabs(ans-oldAns)/ans < 1.0e-4)
      break;
    //    cout << "mMin=" << mMin << " mMax=" << mMax << " fColl=" << ans
    //	 << endl;
    mMin = mMax;
  }
  //  ans /= CRITDENMSOLMPC*c->getOm0hh();
  return ans;
}

// Integrand for collapse fraction.  fcn is massFcn, says which mass function
// to use. 
double setFMCollInt(double mass, Cosmology *c1, double z1, int fcn, 
		  int flag)
{
  static Cosmology *c;
  static double z;
  static int massFcn;
  double result;

  if (flag == 1) {
    c = c1;
    z = z1;
    massFcn = fcn;
    return 0.0;
  }

  if (massFcn == 1)
    result = mass*c->dndlMSheth(z,mass);
  else if (massFcn == 2)
    result = mass*c->dndlMJenkins(z,mass);
  else
    result = mass*c->dndlM(z,mass);
  return result;
}

// Dummy integrand for collapse fraction 
double fMCollInt(double mass)
{
  return setFMCollInt(mass,NULL,0.0,0,0);
}

//////////////////////////////////////////////////////////////////////////
//  Fourier Transform
///////////////////////////////////////////////////////////////////////////
//void FTTable(double xivec[], double rvec[], int nr1, double pvec[], 
//	     double kvec[], int nk1)
///{
//  return;
//}

/*
const double INT_ACCURACY2(1.0e-5);

int nr(100),nk(100);
double *rft,*xift,*xiftsp;

// Function to construct a table of the (3D, isotropic) FT of the
// function input as xivec.  
//   xivec = correlation function values
//   rvec = ordinates of the correlation function (comoving Mpc).
//   nr = number of points in xivec/rvec
//   pvec = vector for storing power spectrum
//   kvec = vector containing wavenumbers at which to evaluate P(k)
//          (comoving Mpc^-1)
//   nk = number of points in kvec/pvec
// Uses spline interpolation on the correlation function to do the 
// evaluation.  Integrates over setFTQRInt.
//
void FTTable(double xivec[], double rvec[], int nr1, double pvec[], 
	     double kvec[], int nk1)
{
  double oldAns;
  double rMin,rMax;
  int ift,br;
  const double natural(1.0e30);
  double npoints(2.0);

  nr=nr1;
  nk=nk1;

  // Set up the spline fit to the correlation function
  rft = dvector(1,nr);
  xift = dvector(1,nr);
  xiftsp = dvector(1,nr);


  for (ift=1;ift<=nr;ift++) {
    rft[ift] = log(rvec[ift]);
    if(xivec[ift]<=0.0){
      cout <<"Non-positive xi in FTTabel"<<endl;
      xift[ift] = log(1.0e-20);  //push to zero     
    }else{
      xift[ift] = log(xivec[ift]);
    }
  }

  spline(rft,xift,nr,natural,natural,xiftsp);
  cout<<"spline OK"<<endl;

  for (ift=1;ift<=nk;ift++) {
    setFTQRInt(0.0,kvec[ift],1);
    rMin = rvec[1];
    rMax = 2.0*PI/kvec[ift]/npoints;
    oldAns = 0.0;
    br=0;
    while (1) {
      // We need to take care that we don't step too far, because 
      // integrand oscillates rapidly 
      if (rMax > rvec[nr]) {
	rMax = rvec[nr]*0.99999;
	br=1;
      }
      oldAns +=  qromb(FTQRInt,rMin,rMax);
      rMin = rMax;
      rMax += 2.0*PI/kvec[ift]/npoints;
      if (br==1) break;
    }
    pvec[ift] = 4.0*PI*oldAns;
  }

  free_dvector(rft,1,nr);
  free_dvector(xift,1,nr);
  free_dvector(xiftsp,1,nr);
  return;
}


// Integrand for making the power spectrum, just a 3D isotropic Fourier
// transform
double setFTQRInt(double r, double kp, int flag)
{
  static double k;
  double xi;
  double ans;
  double lr,lxi;
  if (flag == 1) {
    k = kp;
    return 0.0;
  }
  lr=log(r);
  splint(rft,xift,xiftsp,nr,lr,&lxi);
  xi=exp(lxi);
  ans = r*sin(k*r)/k*xi;
  return ans;
}

// Dummy integrand for making power spectrum 
double FTQRInt(double r)
{
  return setFTQRInt(r,0.0,0);
}
*/

const double INT_ACCURACY2(1.0e-5);

int nr(100),nk(100);
double *rft,*xift,*xiftsp;

// Function to construct a table of the (3D, isotropic) FT of the
// function input as xivec.  
//   xivec = correlation function values
//   rvec = ordinates of the correlation function (comoving Mpc).
//   nr = number of points in xivec/rvec
//   pvec = vector for storing power spectrum
//   kvec = vector containing wavenumbers at which to evaluate P(k)
//          (comoving Mpc^-1)
//   nk = number of points in kvec/pvec
// Uses spline interpolation on the correlation function to do the 
// evaluation.  Integrates over setFTQRInt.
//
void FTTable(double xivec[], double rvec[], int nr1, double pvec[], 
	     double kvec[], int nk1)
{
  double oldAns;
  double rMin,rMax;
  int ift,br;
  const double natural(1.0e30);
  double npoints(2.0);

  nr=nr1;
  nk=nk1;

  // Set up the spline fit to the correlation function
  rft = dvector(1,nr);
  xift = dvector(1,nr);
  xiftsp = dvector(1,nr);


  for (ift=1;ift<=nr;ift++) {
    rft[ift] =(rvec[ift]);
    xift[ift] = (xivec[ift]*rvec[ift]);
  }

  spline(rft,xift,nr,natural,natural,xiftsp);
  cout<<"spline OK"<<endl;

  for (ift=1;ift<=nk;ift++) {
    setFTQRInt(0.0,kvec[ift],1);
    rMin = rvec[1];
    rMax = 2.0*PI/kvec[ift]/npoints;
    oldAns = 0.0;
    br=0;
    while (1) {
      // We need to take care that we don't step too far, because 
      // integrand oscillates rapidly 
      if (rMax > rvec[nr]) {
	rMax = rvec[nr]*0.99999;
	br=1;
      }
      oldAns +=  qromb(FTQRInt,rMin,rMax);
      rMin = rMax;
      rMax += 2.0*PI/kvec[ift]/npoints;
      if (br==1) break;
    }
    pvec[ift] = 4.0*PI*oldAns;
  }

  free_dvector(rft,1,nr);
  free_dvector(xift,1,nr);
  free_dvector(xiftsp,1,nr);
  return;
}


// Integrand for making the power spectrum, just a 3D isotropic Fourier
// transform
double setFTQRInt(double r, double kp, int flag)
{
  static double k;
  double xi;
  double ans;
  double lr,lxi;
  if (flag == 1) {
    k = kp;
    return 0.0;
  }
  lr=r;
  splint(rft,xift,xiftsp,nr,lr,&lxi);
  xi=lxi;
  ans = sin(k*r)/k*xi;
  return ans;
}

// Dummy integrand for making power spectrum 
double FTQRInt(double r)
{
  return setFTQRInt(r,0.0,0);
}


//////////////////////////////////////////////////////////////////////////
// Calculate halo averaged bias
/////////////////////////////////////////////////////////////////////////

//Calculate the bias averaged over mass at a redshift z
// assumes a PS mass function
double TwentyOneCM::meanBiasPS(double z)
{
  double mB(0.0);
  double mN(0.0);
  double mUse,mMin,mMax;
  double v1,v2,dlm;
  int nstep,i;

  //average bias 
  mMin=coolMass(c,z);
  mMax=30.0*mMin;   // Halos should cluster close on cooling mass
  dlm=0.1;
  nstep=int(ceil((log(mMax)-log(mMin))/dlm));

  //integrate weighted bias using trapezium rule
  mUse=mMin;
  v1=nm(mUse,z,c)*biasm(mUse,z,c);
  for(i=0;i<nstep;i++){
    mUse=exp(log(mUse)+dlm);
    v2=nm(mUse,z,c)*biasm(mUse,z,c);
    mB+= 0.5*(v2+v1)*dlm;
    v1=v2;    
    //    cout <<mMin<<"\t"<<mMax<<"\t"<<mUse<<"\t"<<biasm(mUse,z,c)<<endl;

  }
  //integrate number using trapezium rule
  mUse=mMin;
  v1=nm(mUse,z,c);
  for(i=0;i<nstep;i++){
    //mUse+=dlm;
    mUse=exp(log(mUse)+dlm);
    v2=nm(mUse,z,c);
    mN+= 0.5*(v2+v1)*dlm;
    v1=v2;    
  }

  mB/=mN;

  return mB;
}

//Calculate the bias averaged over mass at a redshift z
// assumes a ST mass function
double TwentyOneCM::meanBiasST(double z)
{
  double mB(0.0);
  double mN(0.0);
  double mUse,mMin,mMax;
  double v1,v2,dlm;
  int nstep,i;

  //average bias 
  mMin=coolMass(c,z);
  mMax=30.0*mMin;   // Halos should cluster close on cooling mass
  dlm=0.1;
  nstep=int(ceil((log(mMax)-log(mMin))/dlm));

  //integrate weighted bias using trapezium rule
  mUse=mMin;
  v1=nmST(mUse,z,c)*biasmST(mUse,z,c);
  for(i=0;i<nstep;i++){
    mUse=exp(log(mUse)+dlm);
    v2=nmST(mUse,z,c)*biasmST(mUse,z,c);
    mB+= 0.5*(v2+v1)*dlm;
    v1=v2;    
    //    cout <<mMin<<"\t"<<mMax<<"\t"<<mUse<<"\t"<<biasm(mUse,z,c)<<endl;

  }
  //integrate number using trapezium rule
  mUse=mMin;
  v1=nmST(mUse,z,c);
  for(i=0;i<nstep;i++){
    //mUse+=dlm;
    mUse=exp(log(mUse)+dlm);
    v2=nmST(mUse,z,c);
    mN+= 0.5*(v2+v1)*dlm;
    v1=v2;    
  }

  mB/=mN;

  return mB;
}

//create and use a spline of the Window function
double TwentyOneCM::splineBias(double zuse)
{
  static double *y2;
  static double *xs, *ys;
  static int nmaxs;
  static int icount;
  static int iflag;
  int n,nmax(120);
  double yout;
  double yp1(1.0e30);
  double ypn(1.0e30);
  double xuse;
  char *file="./bias_spline.dat";
  ifstream fin(file);
  ofstream fout;
  
  xuse=zuse;

  //initialise spline on first pass
  if(iflag==0){
    if(icount>0){
      free_dvector(xs,1,nmaxs);
      free_dvector(ys,1,nmaxs);
      free_dvector(y2,1,nmaxs);
    }
    // store information in static variables
    nmaxs=nmax;
    y2=dvector(1,nmax);
    xs=dvector(1,nmax);
    ys=dvector(1,nmax);
    cout <<"initialisaing bias spline"<<endl;
    //if can't initiate from file do calculation
    if(!fin){
      //////////
      double dz(0.1);
      double z;
      //      double zmin(5.0);
      double zmin(ZLOW);
      double zmax(50.0);
      //double zmin(18.0);
      //      double zmax(30.0);
      dz=(zmax-zmin)/double(nmax-1);
      
      cout<<"calculating bias spline"<<endl;
      for(n=1;n<=nmaxs;n++){
	z=zmin+dz*(double)(n-1);
	xs[n]=z;
	if(MASSFCN==1){
	  ys[n]=meanBiasST(z);
	}else{
	  ys[n]=meanBiasPS(z);
	}
      }
      spline(xs,ys,nmaxs,yp1,ypn,y2);
      fout.open(file);
      for(n=1;n<=nmaxs;n++){
	fout << xs[n]<<"\t"<<ys[n]<<"\t"<<y2[n]<<endl;
      }
      fout.close();
    }else{
      for(n=1;n<=nmaxs;n++){
	fin >>xs[n] >>ys[n] >>y2[n];
      }
      fin.close();
    }
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    icount++;
    iflag++;
    cout <<"bias spline initiated"<<endl;
    return yout;
  }
  // apply cubic spline
  if(iflag==1){
    if((xuse<xs[1])||(xuse>xs[nmaxs])) {
      cout <<"x exceeds limits in Bias spline"<<endl;
      return 0.0;
    }
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    return yout;
  }

  if(iflag==2){
    free_dvector(xs,1,nmaxs);
    free_dvector(ys,1,nmaxs);
    free_dvector(y2,1,nmaxs);
    icount=0;
    return 0.0;
  }

  cout << "Error in bias spline"<<endl;
  return 0.0;
}

////////////////////////////////////////////////////////////////////
//  FLUCTUATION POWER SPECTRUM FROM INHOMOGENEOUS XRAY HEATING
////////////////////////////////////////////////////////////////////
// Code to calculate the brightness temperature fluctuations 
// from inhomogeneous Xray heating


//calculate 21cm brightness fluctuations from Xray heating variation
//find delta_T from full integration
double TwentyOneCM::powerSpectrumXray(double z, double k, double store[])
{
  double *betaT, *result;
  double tk,xi,xe,lyaflux;
  double pk,pkfull,ps,psfull;
  double DZ,rFcut,rTcut,rF,rT;
  double tb,gammaA;
  double *gT;

  result=dvector(1,3);
  betaT=dvector(1,5);
  gT=dvector(1,4);

  //calculate thermodynamic quantities
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lyaflux=a->lyaFlux(z);
  tb=(1.0-xi)*tBrightGen(z,tk,xe,lyaflux);
  getBetaTb(z,tk,xi,lyaflux,betaT);

  //calculate delta_T/delta
  gammaA=getGT(z,k,gT);

  //get window function
  //  wk=getWindowXray(z,k);

  //cutoffs and power spectra
  DZ=c->growthFac(z);
  rF=getCutoffFilter(z);
  rT=getCutoffThermal(z,tk);
  rFcut=pow(1.0+k*k*rF*rF,-2.0);
  rTcut=exp(-k*k*rT*rT);

  //I calculate gammaA fully for my main result, but output
  //the uniform heating case for comparison as well

  ps=c->powerSpectrum(k)*DZ*DZ*rFcut*rTcut;
  pk=2.0*ps;
  pk*=betaT[5]*(betaT[1]+betaT[4]*gammaA);

  psfull=k*k*k*2.0*betaT[5]*(betaT[1]+betaT[4]*gT[2])*ps/2.0/PI/PI; 
  pkfull=k*k*k*pk/2.0/PI/PI;

  //output full results for curiosity
  store[1]=gammaA;
  store[2]=2.0*betaT[5]*(betaT[1]+betaT[4]*gT[2])*ps;
  store[3]=pk;
  store[4]=psfull;
  store[5]=pkfull;

  free_dvector(result,1,3);
  free_dvector(betaT,1,5);
  free_dvector(gT,1,4);

  return pkfull;
}



//Calculate the density contribution to the lya fluctuations
double TwentyOneCM::getWindowXray(double z, double k)
{
  double wk(0.0);
  double zmax,zmax2(ZSFR);
  double npoint(2.0);
  double rzmax,zuse;
  double period_max(100.0);
  double b;
  double zcut(ZSFR);
 
  double z2;
  double tk,xe, *result;
  double tol(1.0e-4);
  double zstep;
  int nstep(10);  // setting this larger than 1 causes numerical problems 100
  double lambda(1.0);
  double ldz;

  if(z>ZSFR-1.0) return 0.0; //no x-rays no fluctuations in x-ray heating
  if(k>KCUT) return 0.0; //on small scales uniform heating
  if(z<a->getZReion()) return 0.0;

  result=dvector(1,3);
  a->getTIGMcont(z,result);
  tk=result[1];
  xe=result[3];
  free_dvector(result,1,3);

   //Calculate w(k)
  wk=0.0;
  setWindowXrayK(z,z,k,tk,c,a,this,1);
  setDJxrayDzDE(20.0,z,0.0,c,a,1);
  getDLambdaXrayDzDE(xe,z,a,1);
  //Calculate redshift limit for this line
  zstep=(zmax2-z)/((double)nstep);

  // Sanity check limits of integration: calculate step size
  // in comoving distance, but implement in redshift.
  zmax2=min(z+8.0,ZSFR);
  rzmax=period_max*2.0*PI/k;
  if(rzmax<c->relativeConformalR(z,zmax2)) zcut=c->zOfR(rzmax,z);
  if(fabs(zcut-z)<1.0e-6) zcut=2.0*z;

  // integrate using trapezium rule
  zuse=z;
  z2=c->zOfR(2.0*PI/k/npoint,zuse);
  b=0;
  //break integral into smaller steps
  ldz=exp(log(zmax2/z)/(double)(nstep-1));
  zmax=z*ldz;

  while(1){
    if(z2>zmax){
      if(zmax>=zmax2) b=1;  //finish at edge of star formation
      z2=zmax;
      zmax*=ldz;
    }
    if(z2>zcut){
      //	cout<<"hi"<<endl;
	z2=zcut;
	b=1;
    }
    wk+=qromb(intWindowXrayK,zuse,z2,tol);
    nmarkr++;
    if(b==1) break;
    zuse=z2;
    z2=c->zOfR(2.0*PI/k/npoint,zuse);
  }

  //wk*=c->nh(z)*xe;
  //Calculate total heating rate
  lambda=a->xrayHeating(z,xe);
  wk/=lambda;
  cout << z <<"\t"<< k <<"\t"<<nmarkr<<"\t"<<nmarkt<< "\t" <<wk<<"\t"<<zcut<<endl;
  nmarkt=0;
  nmarkr=0;
  
  return wk;
}

extern double zpsave;  //zpsave is from astrophysics class

//function for integrating window
double setWindowXrayK(double z1, double zuse, double k1,double tk1,Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag)
{
  static Cosmology *c;
  static Astrophysics *a;
  static TwentyOneCM *tocm;
  static double z;
  static double k;
  static double tk;
  double wk,rz,x,oscfun;

  double Emin(RYDBERG);
  double Emax(XrayEmax);
  double dldz(1.0);
  double lEmin,lEmax;
  double tol(1.0e-4);

  lEmin=log(Emin);
  lEmax=log(Emax);
  
  if(iflag==1){
    c=c1;
    a=a1;
    tocm=tocm1;
    z=z1;
    k=k1;
    tk=tk1; //gas kinetic temperature
    return 0.0;
  }

  if(fabs(zuse-z)<1.0e-7){
    // r=0 limit
    oscfun=1.0+tocm->splineBias(zuse);
  }else{
    rz=SPEEDOFLIGHT_CGS*(c->confTime(z)-c->confTime(zuse))/H0_CMSMPC/c->getH();
    x=rz*k;
    oscfun = (1.0+tocm->splineBias(zuse))*sin(x)/x;
    oscfun -= 2.0*((3.0/pow(x,3.0)-1.0/x)*sin(x)-3.0*cos(x)/x/x)/3.0;
  }
  oscfun*=c->growthFac(zuse)/c->growthFac(z);

  //get total X-ray differential flux per atom
  zpsave=zuse;
  dldz=qromb(getDLambdaXrayDlE,lEmin,lEmax,tol);

  wk=dldz*oscfun;
 
  nmarkt++;

  return wk;
  // return dldz; //if return dldz then should get wk=1
}

//
double intWindowXrayK(double zuse)
{
  return setWindowXrayK(1.0,zuse,0.0,0.0,NULL,NULL,NULL,0);
}

double TwentyOneCM::powerSpectrumFull(double z, double k, double store[], int iflag)
{
  double *betaT, *result;
  double tk,xi,xe,lyaflux;
  double pkfull,ps;
  double DZ,rFcut,rTcut,rF,rT;
  double tb,gammaA(0.0),wklya(0.0),wkxray(0.0);
  double *gT;
  double pkiso,pk0,pk2,pk4;
  double k3p;
  double ge(0.0),gX;

  result=dvector(1,3);
  betaT=dvector(1,5);
  gT=dvector(1,4);

  //calculate thermodynamic quantities
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lyaflux=a->lyaFlux(z);
  tb=(1.0-xi)*tBrightGen(z,tk,xe,lyaflux);
  getBetaTb(z,tk,xi,lyaflux,betaT);

  //calculate delta_T/delta and W_\alpha(k)
  if(iflag==0){  //adiabatic cooling no lya fluct
    getGT(z,k,gT);
    gammaA=gT[3];
    wklya=0.0;
    wkxray=0.0;
  }else   if(iflag==1){  //full case
    getGT(z,k,gT);
    gammaA=gT[1]; 
    wklya=getWindowLya(z,k);
    wkxray=getWindowXray(z,k);
  }else if(iflag==2){  //X-ray fluct no lya
    getGT(z,k,gT);
    gammaA=gT[1]; 
    wklya=0.0;
    wkxray=getWindowXray(z,k);
  }else if(iflag==3){  //lya fluct uniform heating
    getGT(z,k,gT);
    gammaA=gT[2]; 
    wklya=getWindowLya(z,k);
    wkxray=0.0;
  }else if(iflag==4){  //no lya uniform heating
    getGT(z,k,gT);
    gammaA=gT[2]; 
    wklya=0.0;
    wkxray=0.0;
  }else if(iflag==5){  //no lya - lambda approx
    //    getGT(z,k,gT);
    wkxray=getWindowXray(z,k);
    gammaA=wkxray; 
    wklya=0.0;
  }else if(iflag==6){  //lya+xray - lambda approx
    wklya=getWindowLya(z,k);
    wkxray=getWindowXray(z,k);
    gammaA=wkxray;
  }else if(iflag==7){  //lya+adiabatic cooling
    getGT(z,k,gT);
    gammaA=gT[3]; 
    wklya=getWindowLya(z,k);
    wkxray=0.0;
  }else   if(iflag==8){  //full case +gX
    getGT(z,k,gT);
    gammaA=gT[1];
    ge=gT[4];
    wklya=getWindowLya(z,k);
    wkxray=getWindowXray(z,k);
  }

  gX=-xe/(1.0-xe)*ge;

  //cutoffs and power spectra
  DZ=c->growthFac(z);
  rF=getCutoffFilter(z);
  rT=getCutoffThermal(z,tk);
  rFcut=pow(1.0+k*k*rF*rF,-2.0);
  rTcut=exp(-k*k*rT*rT);

  //calculate 3 power spectra
  ps=c->powerSpectrum(k)*DZ*DZ*rFcut*rTcut;

  pkiso=betaT[1]+betaT[2]*gX+betaT[3]*wklya+betaT[4]*gammaA;
  pk0=pkiso*pkiso;
  pk2=2.0*pkiso*betaT[5];
  pk4=betaT[5]*betaT[5];

  pk0*=ps;
  pk2*=ps;
  pk4*=ps;
  pkiso=pk0+pk2/3.0+pk4/5.0;

  k3p=k*k*k/2.0/PI/PI;
  pkfull=sqrt(pkiso*k3p);

  //output full results for curiosity
  store[1]=pk0;
  store[2]=pk2;
  store[3]=pk4;
  store[4]=pkiso;
  store[5]=gammaA;
  store[6]=wklya;
  store[7]=wkxray;
  store[8]=gX;

  free_dvector(result,1,3);
  free_dvector(betaT,1,5);
  free_dvector(gT,1,4);

  return pkfull;
}

// Calculate nz power spectrum at redshifts zin[nz] at wavenumber k and 
// return in store[7][nz]

// iflag option:
// 0: adiabatic cooling + no lya
// 1: xray fluc + lya 
// 2: xray fluc + no lya
// 3: uniform + lya
// 4: uniform + no lya
// 5: lambda approx + no lya
// 6: lambda approx + lya
// 7: adiabatic + lya
// 8: gT=2/3 + lya
// 9: xrayfluc + lya +gX

void TwentyOneCM::powerSpectrumFullN(double zin[], double k, double **store,int nz, int iflag)
{
  double *betaT, *result;
  double tk,xi,xe,lyaflux;
  double pkfull,ps;
  double DZ,rFcut,rTcut,rF,rT;
  double tb,gammaA(0.0),wklya(0.0),wkxray(0.0);
  double **gT;
  double pkiso,pk0,pk2,pk4;
  double k3p;
  int i;
  double z;
  double ge(0.0),gX;

  result=dvector(1,3);
  betaT=dvector(1,5);
  gT=dmatrix(1,4,1,nz);

  //appropraite initialisation
  if(iflag!=5 && iflag!=6 && iflag!=8) getGTN(zin,k,gT,nz);

  //loop over redshifts
  for(i=1;i<=nz;i++){
    z=zin[i];
 
    //calculate delta_T/delta and W_\alpha(k)
    if(iflag==0){  //adiabatic cooling no lya fluct
      gammaA=gT[3][i];
      wklya=0.0;
      wkxray=0.0;
    }else   if(iflag==1){  //full case
      gammaA=gT[1][i]; 
      wklya=getWindowLya(z,k);
      wkxray=splineWkXray(z,k,1);
    }else if(iflag==2){  //X-ray fluct no lya
      gammaA=gT[1][i]; 
      wklya=0.0;
      wkxray=splineWkXray(z,k,1);
    }else if(iflag==3){  //lya fluct uniform heating
      gammaA=gT[2][i]; 
      wklya=getWindowLya(z,k);
      wkxray=0.0;
    }else if(iflag==4){  //no lya uniform heating
      gammaA=gT[2][i]; 
      wklya=0.0;
      wkxray=0.0;
    }else if(iflag==5){  //no lya - lambda approx
      wkxray=getWindowXray(z,k);
      gammaA=wkxray; 
      wklya=0.0;
    }else if(iflag==6){  //lya+xray - lambda approx
      wklya=getWindowLya(z,k);
      wkxray=getWindowXray(z,k);
      gammaA=wkxray;
    }else if(iflag==7){  //lya+adiabatic cooling
      gammaA=gT[3][i]; 
      wklya=getWindowLya(z,k);
      wkxray=0.0;
    }else if(iflag==8){  //lya+g=2/3
      gammaA=2.0/3.0; 
      wklya=getWindowLya(z,k);
      wkxray=0.0;
    }else if(iflag==9){  //full case + gX
      gammaA=gT[1][i]; 
      ge=gT[4][i];
      wklya=getWindowLya(z,k);
      wkxray=splineWkXray(z,k,1);
    }

    //calculate thermodynamic quantities
    a->getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    lyaflux=a->lyaFlux(z);
    tb=(1.0-xi)*tBrightGen(z,tk,xe,lyaflux);
    getBetaTb(z,tk,xi,lyaflux,betaT);
    
    gX=-xe/(1.0-xe)*ge;

    //cutoffs and power spectra
    DZ=c->growthFac(z);
    rF=getCutoffFilter(z);
    rT=getCutoffThermal(z,tk);
    rFcut=pow(1.0+k*k*rF*rF,-2.0);
    rTcut=exp(-k*k*rT*rT);
    
    //calculate 3 power spectra
    ps=c->powerSpectrum(k)*DZ*DZ*rFcut*rTcut;
    
    pkiso=betaT[1]+betaT[2]*gX+betaT[3]*wklya+betaT[4]*gammaA;
    pk0=pkiso*pkiso;
    pk2=2.0*pkiso*betaT[5];
    pk4=betaT[5]*betaT[5];
    
    pk0*=ps;
    pk2*=ps;
    pk4*=ps;
    pkiso=pk0+pk2/3.0+pk4/5.0;
    
    k3p=k*k*k/2.0/PI/PI;
    pkfull=sqrt(pkiso*k3p);

    //output full results for curiosity
    store[1][i]=pk0;
    store[2][i]=pk2;
    store[3][i]=pk4;
    store[4][i]=pkiso;
    store[5][i]=gammaA;
    store[6][i]=wklya;
    store[7][i]=wkxray;
    store[8][i]=gX;
    
  }
  free_dvector(result,1,3);
  free_dvector(betaT,1,5);
  free_dmatrix(gT,1,4,1,nz);
  
}

////////////////////////////////////////////////////////////////////
//  Temperature power spectrum
////////////////////////////////////////////////////////////////////


double TwentyOneCM::powerSpectrumTemp(double z, double k, double store[])
{
  double *result;
  double tk,xi,xe;
  double ps;
  double DZ,rFcut,rTcut,rF,rT;
  double wkxray(0.0);
  double *gT;
  double pkA,pkU,pkT;
  double k3p;

  result=dvector(1,3);
  gT=dvector(1,4);

  //calculate thermodynamic quantities
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];

  //calculate delta_T/delta 
  getGT(z,k,gT);
  wkxray=getWindowXray(z,k);

  //cutoffs and power spectra
  DZ=c->growthFac(z);
  rF=getCutoffFilter(z);
  rT=getCutoffThermal(z,tk);
  rFcut=pow(1.0+k*k*rF*rF,-2.0);
  rTcut=exp(-k*k*rT*rT);

  //calculate 3 power spectra
  ps=c->powerSpectrum(k)*DZ*DZ*rFcut*rTcut;

  pkA=gT[3]*gT[3]*tk*tk*ps;
  pkU=gT[2]*gT[2]*tk*tk*ps;
  pkT=gT[1]*gT[1]*tk*tk*ps;

  k3p=k*k*k/2.0/PI/PI;

  //output full results for curiosity
  store[1]=pkT*k3p;
  store[2]=pkU*k3p;
  store[3]=pkA*k3p;
  store[4]=gT[1];
  store[5]=gT[2];
  store[6]=gT[3];
  store[7]=wkxray;
  store[8]=tk;

  free_dvector(result,1,3);
  free_dvector(gT,1,4);

  return pkT*k3p;
}

void TwentyOneCM::powerSpectrumTempN(double zin[], double k, double **store, int nz)
{
  double *result;
  double tk,xi,xe;
  double ps;
  double DZ,rFcut,rTcut,rF,rT;
  double wkxray(0.0);
  double **gT;
  double pkA,pkU,pkT;
  double k3p;
  int i;
  double z;

  result=dvector(1,3);
  gT=dmatrix(1,4,1,nz);

  //calculate delta_T/delta 
  getGTN(zin,k,gT,nz);

  for(i=1;i<=nz;i++){
    z=zin[i];
    //calculate thermodynamic quantities
    a->getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    
    wkxray=getWindowXray(z,k);
    
    //cutoffs and power spectra
    DZ=c->growthFac(z);
    rF=getCutoffFilter(z);
    rT=getCutoffThermal(z,tk);
    rFcut=pow(1.0+k*k*rF*rF,-2.0);
    rTcut=exp(-k*k*rT*rT);
    
    //calculate 3 power spectra
    ps=c->powerSpectrum(k)*DZ*DZ*rFcut*rTcut;
    
    pkA=gT[3][i]*gT[3][i]*tk*tk*ps;
    pkU=gT[2][i]*gT[2][i]*tk*tk*ps;
    pkT=gT[1][i]*gT[1][i]*tk*tk*ps;
    
    k3p=k*k*k/2.0/PI/PI;
    
    //output full results for curiosity
    store[1][i]=pkT*k3p;
    store[2][i]=pkU*k3p;
    store[3][i]=pkA*k3p;
    store[4][i]=gT[1][i];
    store[5][i]=gT[2][i];
    store[6][i]=gT[3][i];
    store[7][i]=wkxray;
    store[8][i]=tk;
    
  }

  free_dvector(result,1,3);
  free_dmatrix(gT,1,4,1,nz);

}


////////////////////////////////////////////////////////////////////
//  Evolution of Temperature perturbations
////////////////////////////////////////////////////////////////////

double TwentyOneCM::QXray(double z)
{
  double Q,ntot,CT,dzdt,heat;
  double tk,xi,xe,*result;
  double zz;

  result=dvector(1,3);
  zz=z+1.0;
  a->getTIGMcont(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  ntot=c->nh(z)*(1.0+xe+FHE);  //n=n_H+n_e+n_He
  CT=1.5*BOLTZK*ntot;  //specific heat capacity in ergs K^{-1} cm^{-3}
  dzdt=zz*c->hubbleZ(z)*UNH*c->getH();
  
  //heat=xrayHeating(z); //xrays
  heat=a->lumTotStarBurst(a->globalSFR(z))/MPC/MPC/MPC*a->fracPrimaryElectron(xe,1); //xrays
  heat/=dzdt*CT; 
  
  Q=heat/tk;  //Coefficient on RHS of delta_T evolution equation
  free_dvector(result,1,3);

  return Q;
}

double TwentyOneCM::QCompton(double z)
{
  double Q,dzdt,heat;
  double tk,xi,xe,*result;
  double zz,tgammainv,Tcmb,ucmb;
  

  result=dvector(1,3);
  zz=z+1.0;
  a->getTIGMcont(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  dzdt=zz*c->hubbleZ(z)*UNH*c->getH();
  
  Tcmb=TCMB*(1.0+z);
  ucmb=RADIATIONA*Tcmb*Tcmb*Tcmb*Tcmb;
  tgammainv=8.0*ucmb*SIGMATHOMSON/3.0/ELECTRONMASS/SPEEDOFLIGHT_CGS;
  
  heat=xe/(1.0+xe+FHE);
  heat*=tgammainv;
  heat*=Tcmb/tk;
  
  Q=heat/dzdt;  //Coefficient on RHS of delta_T evolution equation
  free_dvector(result,1,3);

  return Q;
}

double TwentyOneCM::QIon(double z)
{
  double Q,dzdt,ion;
  double tk,xi,xe,*result;
  double zz;

  result=dvector(1,3);
  zz=z+1.0;
  a->getTIGMcont(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  dzdt=zz*c->hubbleZ(z)*UNH*c->getH();
  
  //heat=xrayHeating(z); //xrays
  ion=a->lumTotStarBurst(a->globalSFR(z))/MPC/MPC/MPC*a->fracPrimaryElectron(xe,2); //xrays
  ion/=RYDBERG_CGS;
  ion/=c->nh(z);  // ionization rate per baryon
  ion*=(1-xe)/xe;
  ion/=dzdt; 
  
  Q=ion;  //Coefficient on RHS of delta_x evolution equation
  free_dvector(result,1,3);

  return Q;
}

double TwentyOneCM::QRec(double z)
{
  double Q,dzdt,rec;
  double tk,xi,xe,*result;
  double zz;
  double alphaA(4.2e-13);  //case-A recombination: cm^3 s^-1 (T=10^4K)

  result=dvector(1,3);
  zz=z+1.0;
  a->getTIGMcont(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  dzdt=zz*c->hubbleZ(z)*UNH*c->getH();
  
  rec=alphaA*CLUMPING*xe*c->nh(z);
  rec/=dzdt; 
  
  Q=rec;
  free_dvector(result,1,3);

  return Q;
}


//evolve g_T=delta_T/delta
double TwentyOneCM::getGTHistory(double zin, double k)
{
  double zstart(1000.0),zend;

  double *y;
  double *dy;
  int nvar(4);
  int nok,nbad;
  double hstep(1.0);
  double hmin(0.001);
  double eps(1e-4);
  double gT,QC,QX,QI,QR;
  ofstream fout;
  char *file="gT.dat";
  double tk,xi,xe,*result;
  double gXH;
  
  y=dvector(1,nvar);
  dy=dvector(1,nvar);
  result=dvector(1,4);

  setDerivsGT(zin,k,y,dy,c,this,1);
  splineWkXray(10.0,k,0);

  y[1]=0.0;  //full
  y[2]=0.0;  //uniform
  y[3]=0.0;  //compton
  y[4]=0.0;  //x_e

  fout.open(file);
  zend=zin;
  zend=zstart-1.0;
  while(zend>5.0){
    odeint(y,nvar,zstart,zend,eps,hstep,hmin,&nok,&nbad,derivsGT,rkqs);

    a->getTIGM(zend,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];

    gT=y[1];
    QC=QCompton(zend);
    QX=QXray(zend);
    QI=QIon(zend);
    QR=QRec(zend);
    gXH=y[4]*-xe/(1.0-xe);
    cout<<zend<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<y[3]<<"\t"<<y[4]<<"\t"<<gXH<<"\t"<<QC<<"\t"<<QX<<"\t"<<QI<<"\t"<<QR<<endl;
    fout<<zend<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<y[3]<<"\t"<<y[4]<<"\t"<<gXH<<"\t"<<QC<<"\t"<<QX<<"\t"<<QI<<"\t"<<QR<<endl;

    zstart=zend;
    zend-=1.0;
  }
  gT=y[1];
  fout.close();

  free_dvector(y,1,nvar);  
  free_dvector(dy,1,nvar);
  free_dvector(result,1,4);

  return gT;
}

//evolve g_T=delta_T/delta
double TwentyOneCM::getGT(double zin, double k, double gT[])
{
  double zstart(1000.0),zend;

  double *y;
  double *dy;
  int nvar(4);
  int nok,nbad;
  double hstep(1.0);
  double hmin(0.001);
  double eps(1e-4);
  double gTFull;
  double zmin(10.0);

  y=dvector(1,nvar);
  dy=dvector(1,nvar);

  setDerivsGT(0.0,k,y,dy,c,this,1);

  zmin=0.9*zin;
  splineWkXray(zmin,k,0);

  y[1]=0.0;
  y[2]=0.0;
  y[3]=0.0;
  y[4]=0.0;

  zend=zin;
  odeint(y,nvar,zstart,zend,eps,hstep,hmin,&nok,&nbad,derivsGT,rkqs);
  
  gTFull=y[1];
  
  gT[1]=y[1];
  gT[2]=y[2];
  gT[3]=y[3];  
  gT[4]=y[4];


  free_dvector(y,1,nvar);  
  free_dvector(dy,1,nvar);

  return gTFull;
}

//evolve g_T=delta_T/delta
//provide vector zin[nz] with z values to obtain gT. 
//Order from low z to high z 
void TwentyOneCM::getGTN(double zin[], double k, double **gT, int nz)
{
  double zstart(1000.0),zend;

  double *y;
  double *dy;
  int nvar(4);
  int nok,nbad;
  double hstep(1.0);
  double hmin(0.001);
  double eps(1e-4);
  double gTFull;
  int i;
  double zmin(10.0);

  y=dvector(1,nvar);
  dy=dvector(1,nvar);

  setDerivsGT(0.0,k,y,dy,c,this,1);

  zmin=0.9*fmin(zin[1],zin[nz]);

  splineWkXray(zmin,k,0);

  y[1]=0.0;
  y[2]=0.0;
  y[3]=0.0;
  y[4]=0.0;
  ///////////////////
  
  //scale dependent initial conditions
  /*
  zstart=100.0;
  if(k<1.0e-2){
    y[1]=0.68;
    y[2]=0.68;
    y[3]=0.68;
  }else {
    y[1]=0.29;
    y[2]=0.29;
    y[3]=0.29;
  }    
  */ 

  /////////////////

  //  zstart=1000.0;
  for(i=nz;i>=1;i--){
    zend=zin[i];
    odeint(y,nvar,zstart,zend,eps,hstep,hmin,&nok,&nbad,derivsGT,rkqs);
    
    gTFull=y[1];
    
    gT[1][i]=y[1];
    gT[2][i]=y[2];
    gT[3][i]=y[3];
    gT[4][i]=y[4];    
    zstart=zend;
  }

  free_dvector(y,1,nvar);  
  free_dvector(dy,1,nvar);

}

// function to calculate derivatives for globalHistory
void setDerivsGT(double z, double k1, double y[],double dy[],Cosmology *c1,TwentyOneCM *tocm1,int iflag)
{
  static Cosmology *c;
  static TwentyOneCM *tocm;
  double gT,dgT,QC,QX;
  double gX,dgX,QI,QR;
  double wk(0.0);
  static double k;

  if(iflag==1){
    c=c1;
    tocm=tocm1;
    k=k1;
  }else{
    //    wk=tocm->getWindowXray(z,k);
    wk=tocm->splineWkXray(z,k,1);
    QC=tocm->QCompton(z);
    QX=tocm->QXray(z);
    QI=tocm->QIon(z);
    QR=tocm->QRec(z);
    //Full calculation
    gT=y[1];
    dgT=(gT-2.0/3.0)/(1.0+z);
    dgT+=QC*gT;
    dgT-=QX*(wk-gT);
    dy[1]=dgT;
    //uniform heating
    gT=y[2];
    dgT=(gT-2.0/3.0)/(1.0+z);
    dgT+=QC*gT;
    dgT-=QX*(-gT);
    dy[2]=dgT;
    //Compton only
    gT=y[3];
    dgT=(gT-2.0/3.0)/(1.0+z);
    dgT+=QC*gT;
    dy[3]=dgT;

    //x_e fluctuations 
    gX=y[4];
    dgX=gX/(1.0+z);
    dgX-=QI*(wk-gX);
    dgX+=QR*(1.0+gX);
    dy[4]=dgX;

  }

}

//dummy function to get derivative for odeint for globalHistory
void derivsGT(double z, double y[], double dy[])
{
  setDerivsGT(z,0.0,y,dy,NULL,NULL,0);
}



//create and use a spline of the xray Window function
//
// On initialisation zuse specifies lowest z value of spline
//
double TwentyOneCM::splineWkXray(double zuse, double k, int useflag)
{
  static double *y2;
  static double *xs, *ys;
  static int nmaxs;
  static int icount;
  static int iflag;
  int n,nmax(35);
  nmax=25;

  double yout;
  double yp1(1.0e30);
  double ypn(1.0e30);
  double xuse;
  
  xuse=zuse;
  if(xuse>ZSFR) return 0.0;

  //initialise spline on first pass or when useflag==0
  if(useflag==0){
    iflag=0;
    if(icount>0){
      //free_dvector(xs,1,nmaxs);
      //free_dvector(ys,1,nmaxs);
      //free_dvector(y2,1,nmaxs);
    }
    icount=0;
    // store information in static variables
    nmaxs=nmax;
    y2=dvector(1,nmax);
    xs=dvector(1,nmax);
    ys=dvector(1,nmax);
    cout <<"initialisaing xray spline"<<endl;
    //////////
    double dz(0.1);
    double z;
    double zmin(8.0);
    double zmax(ZSFR);

    zmin=zuse;
    dz=(zmax-zmin)/double(nmax-1);
    
    cout<<"calculating xray spline"<<endl;
    for(n=1;n<=nmaxs;n++){
      z=zmin+dz*(double)(n-1);
      xs[n]=z;
      ys[n]=getWindowXray(z,k);
    }
    spline(xs,ys,nmaxs,yp1,ypn,y2);
    
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    icount++;
    iflag++;
    cout <<"xray spline initiated"<<endl;
    return yout;
  }
  // apply cubic spline
  if(iflag==1){
    if((xuse<xs[1])||(xuse>xs[nmaxs])) {
      cout <<"x exceeds limits in xray spline"<<endl;
      cout <<xuse<<endl;
      return 0.0;
    }
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    return yout;
  }

  cout << "Error in xray spline"<<endl;
  return 0.0;
}
//////////////////////////////////////////////////////////////////////
// Poisson fluctuation code - T
//////////////////////////////////////////////////////////////////////

//WORKING UNDER THE delta_T=delta_Lambda approximation

//Calculate the Poisson contribution to the T fluctuations
//because I have to calculate the correlation function and then 
//fourier transform it, its not worth getting any single k value
//I just do everything at once and store the correlation function
//and the power spectrum in files.

double TwentyOneCM::poissonPowerXray(double z)
{
  double pC;
  int  nr(20),nk,nuse(20);
  double l,lmin,lmax,dl;
  double k,kmin,kmax,dk;
  double *lvec,*pcvec;
  double *kvec,*ftvec;
  double tk,xi,xe,lyaflux,*betaT,*result;
  double rT,ps,pscut,psfull,rTcut;
  double smallfac(8.0);
  int i;
  char *file;
  ofstream fout;

  cout<<c->relativeConformalR(z,ZSFR)<<endl;

  //Thermodynamics
  result=dvector(1,3);
  betaT=dvector(1,5);
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lyaflux=a->lyaFlux(z);
  getBetaTb(z,tk,xe,lyaflux,betaT);

  //Set limits of correlation function
  lmax=2.0*c->relativeConformalR(z,ZSFR);  //NEED TO MODIFY
  lmin=1.0e-4;
  dl=exp(log(lmax/lmin)/(double)(nuse));
  cout<<lmin<<"\t"<<lmax<<"\t"<<dl<<"\t"<<nr<<endl;
  //find number of data points taken
  i=0;
  l=lmin;
  while(l<lmax){
    i++;
    if(l>lmax/4.0){ //when near end tighten up steps
      l*=1.0+(dl-1.0)/smallfac;
    }else{
      l*=dl;
    }
  }
  nr=i;
  cout<<lmin<<"\t"<<lmax<<"\t"<<dl<<"\t"<<nr<<endl;
  lvec=dvector(1,nr);
  pcvec=dvector(1,nr);

  //Calculate correlation function
  cout<<"Calculating Xray correlation function"<<endl;
  file="poisson_pc_xray.dat";
  fout.open(file);
  i=0;
  l=lmin;
  while(l<=lmax){
    i++;
    pC=poissonCorrelationXray(l,z);
    lvec[i]=l;
    pcvec[i]=pC;

    fout <<l<<"\t"<<pC<<endl;
    cout<<z<<"\t"<<l<<"\t"<<pC<<endl;

    if(l>lmax/4.0){ //when near end tighten up steps
      l*=1.0+(dl-1.0)/smallfac;
    }else{
      l*=dl;
    }
  }
  fout.close();

  //Now perform Fourier Transform
  nk=2*nuse;
  kvec=dvector(1,nk);
  ftvec=dvector(1,nk);

  kmax=1.0e3;;
  kmin=1.0e-3;
  dk=exp(log(kmax/kmin)/(double)(nk-1));

  //set kvector
  k=kmin;
  for(i=1;i<=nk;i++){
    kvec[i]=k;
    k*=dk;
  }
  cout <<"doing FT"<<endl;
  FTTable(pcvec,lvec,nr,ftvec,kvec,nk);
  cout<<"FT done"<<endl;

  file="poisson_ft_xray.dat";
  fout.open(file);
  for(i=1;i<=nk;i++){
    fout<<kvec[i]<<"\t"<<ftvec[i]<<endl;
   }
  fout.close();

  cout<<"calculating Xray power spectrum"<<endl;
  //calculate the cutoff scales
  rT=8.3*sqrt(tk/100.0)/sqrt((1.0+z)/7)*1e-3;
  file="poisson_tb_xray.dat";
  fout.open(file);
  fout<<"#"<<rT<<endl;
  for(i=1;i<=nk;i++){
    k=kvec[i];
    ps=ftvec[i];
    rTcut=exp(-k*k*rT*rT);
    pscut=ps*rTcut;
    pscut*=betaT[4]*betaT[4];
    psfull=(k*k*k*pscut/2.0/PI/PI);
    cout << k <<"\t"<<ps<<"\t"<<pscut<<"\t"<<psfull<<endl;
    fout << k <<"\t"<<ps<<"\t"<<pscut<<"\t"<<psfull<<endl;   
  }
  fout.close();
  cout<<"Xray power spectrum done"<<endl;

  free_dvector(result,1,3);
  free_dvector(betaT,1,5);
  free_dvector(lvec,1,nr);
  free_dvector(pcvec,1,nr);
  free_dvector(kvec,1,nk);
  free_dvector(ftvec,1,nk);

  return rT;
}



double TwentyOneCM::weightPXray(double z,double zp, double tk, double xi)
{
  double wP;
  double dzdr;

  if(z>ZSFR) return 0.0;

  dzdr=c->hubbleZ(zp)*H0_CMSMPC*c->getH()/SPEEDOFLIGHT_CGS;
  wP=getDLambdaXrayDz(zp)*dzdr;
  wP/=4.0*PI*a->splineFColl(zp,c)* CRITDENMSOLMPC*c->getOm0hh();
  return wP;
}

double setPoissonCorrelationKernelXray(double rA, double theta, double z1,double l1, double tk1, double xi1, Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag)
{
  double pC,rB;
  double zpA,zpB;
  static Cosmology *c;
  static Astrophysics *a;
  static TwentyOneCM *tocm;
  static double z;
  static double l;
  static double tk;
  static double xi;
  if(iflag==1){
    c=c1;
    a=a1;
    z=z1;
    l=l1;
    tk=tk1;
    xi=xi1;
    tocm=tocm1;
    return 0.0;
  }

  //impose the boundary limit 
  if(fabs(rA)<1.0e-6) {
    cout <<"using min bound"<<endl;
   return 0.0; 
  }
  if(rA*cos(theta)> l/2.0+1.0e-8)  {
    cout << "rA*cos(theta)>l/2" <<endl;
    //return 0.0;
  }
 
  zpA=c->zOfR(rA,z);
  rB=sqrt(l*l+rA*rA-2.0*l*rA*cos(theta));
  zpB=c->zOfR(rB,z);

  // avoid problems if out of spline reach, at these distances one of
  // P(zA) or P(zB) will be zero as beyond max distance contributing to Xray.
  if(zpA>ZSFR || zpB>ZSFR){
    cout<<"beyond starformation epoch"<<endl;
    return 0.0;
  }


  pC=rA*rA*sin(theta);
  pC*=tocm->weightPXray(z,zpA,tk,xi)*tocm->weightPXray(z,zpB,tk,xi)/rA/rA/rB/rB;
  pC*=a->splineFColl(zpB,c)/a->splineFColl(zpA,c);
  pC*=tocm->splineMFColl(zpA);

  return pC;
}

//dummy function for 2D integration
double intPoissonXray(double rA, double theta)
{
  return setPoissonCorrelationKernelXray(rA,theta,0.0,0.0,0.0,0.0,NULL,NULL,NULL,0);
}

double TwentyOneCM::poissonCorrelationXray(double l, double z)
{
  double pC(0.0),rmax;
  double tol(1.0e-3);
  double lambda;
  double lrmin,lrmax;
  double rxray;
  double tk,xi,xe,*result;
  lsave=l;

  //Basic thermodynamic quantities
  result=dvector(1,3);
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lambda=a->xrayHeating(z,xe);

  setPoissonCorrelationKernelXray(1.0,1.0,z,l,tk,xe,c,a,this,1);
  setDJxrayDzDE(20.0,z,0.0,c,a,1);
  getDLambdaXrayDzDE(xe,z,a,1);
  //limit integral to twice distance to when heating switches on
  rxray=c->relativeConformalR(z,ZSFR);

  //NEED TO CONSIDER NEXT STEP MORE CAREFULLY
  rmax=fmin(10.0*l,rxray*0.99999);
  lrmin=log(1.0e-5);
  //limit integration to region where can get contribution
  if(l>rxray)  lrmin=log(l-rxray+1.0e-5);
  lrmax=log(rmax);
  
  //Calculate correlation function
  pC=qromb(intPoissonRXray,lrmin,log(rmax),tol);
  pC*=2.0*PI;
  pC*=2.0/lambda/lambda;

  cout <<lsave<<"\t"<<"\t"<<nmarkr<<"\t"<<nmarkt<<"\t"<<pC<<endl;
  nmarkr=0;
  nmarkt=0;

  free_dvector(result,1,3);

  return pC;
}

//integral over log r for Poisson fluctuations
double intPoissonRXray(double lrA)
{
  double iPR;
  double tol(1.0e-3);
  double thetamin(0.0);
  double thetamax(PI);
  double rA;
  double rxray;
  
  rA=exp(lrA);
 //rxray=2.0*max distance Xray can contribute
  rxray=2.0*2000.0;  //NEED TO MODIFY!!

  if(rA>rxray) return 0.0;

  //restrict geometry
  if(rA>lsave/2.0){ //restrict to rA<rB
    thetamin=acos(lsave/2.0/rA);
  }
  if(rA+lsave>rxray){ //restrict to overlap region
    thetamax=rA*rA+lsave*lsave-rxray*rxray;
    thetamax/=2.0*rA*lsave;
    if(fabs(thetamax)<=1.0) {
      thetamax=acos(thetamax);
    }else{
      return 0.0;
    }
    if(thetamax<=thetamin) return 0.0;
  }

  rAsave=rA;
  iPR=qromb(intPoissonThetaXray,thetamin,thetamax,tol);
  iPR*=rA;

  // cout <<lsave<<"\t"<<rA<<"\t"<<thetamark<<"\t"<<nmarkr<<"\t"<<nmarkt<<"\t"<<iPR<<endl;
  //nmarkt=0;
  nmarkr++;

  return iPR;
}

double intPoissonThetaXray(double theta)
{
  double iPT;
  nmarkt++;
  iPT= intPoissonXray(rAsave,theta);
  return iPT;
}

//////////////////////////////////////////////////////////////////////
// Radio fluctuation code
//////////////////////////////////////////////////////////////////////
//Calculate the density contribution to the Radio fluctuations
double TwentyOneCM::getWindowRadio(double z, double k)
{
  double wk(0.0);
  double zmax,zmax2(ZSFR);
  double npoint(2.0);
  double rzmax,zuse;
  double period_max(100.0);
  double b;
  double zcut(ZSFR);
 
  double z2;
  double tol(1.0e-4);
  double zstep;
  int nstep(10);  // setting this larger than 1 causes numerical problems 100
  double radio(1.0);
  double ldz;

  if(z>ZSFR) return 0.0; //no x-rays no fluctuations in x-ray heating

   //Calculate w(k)
  wk=0.0;
  setWindowRadioK(z,z,k,c,a,this,1);
  setDJradioDz(nu21cm,z,0.0,c,a,1);

  //Calculate redshift limit for this line
  zstep=(zmax2-z)/((double)nstep);

  // Sanity check limits of integration: calculate step size
  // in comoving distance, but implement in redshift.
  zmax2=min(z+8.0,ZSFR);
  rzmax=period_max*2.0*PI/k;
  if(rzmax<c->relativeConformalR(z,zmax2)) zcut=c->zOfR(rzmax,z);
  if(fabs(zcut-z)<1.0e-6) zcut=2.0*z;

  // integrate using trapezium rule
  zuse=z;
  z2=c->zOfR(2.0*PI/k/npoint,zuse);
  b=0;
  //break integral into smaller steps
  ldz=exp(log(zmax2/z)/(double)(nstep-1));
  zmax=z*ldz;

  while(1){
    if(z2>zmax){
      if(zmax>=zmax2) b=1;  //finish at edge of star formation
      z2=zmax;
      zmax*=ldz;
    }
    if(z2>zcut){
      //	cout<<"hi"<<endl;
	z2=zcut;
	b=1;
    }
    wk+=qromb(intWindowRadioK,zuse,z2,tol);
    nmarkr++;
    if(b==1) break;
    zuse=z2;
    z2=c->zOfR(2.0*PI/k/npoint,zuse);
  }

  //wk*=c->nh(z)*xe;
  //Calculate total heating rate
  radio=a->radioFlux(nu21cm,z);
  wk/=radio;
  cout << z <<"\t"<< k <<"\t"<<nmarkr<<"\t"<<nmarkt<< "\t" <<wk<<"\t"<<zcut<<endl;
  nmarkt=0;
  nmarkr=0;
  
  return wk;
}

//function for integrating window
double setWindowRadioK(double z1, double zuse, double k1,Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag)
{
  static Cosmology *c;
  static Astrophysics *a;
  static TwentyOneCM *tocm;
  static double z;
  static double k;
  double wk,rz,x,oscfun;

  double dJdz(1.0);

  if(iflag==1){
    c=c1;
    a=a1;
    tocm=tocm1;
    z=z1;
    k=k1;
    return 0.0;
  }

  if(fabs(zuse-z)<1.0e-7){
    // r=0 limit
    oscfun=1.0+tocm->splineBias(zuse);
  }else{
    rz=SPEEDOFLIGHT_CGS*(c->confTime(z)-c->confTime(zuse))/H0_CMSMPC/c->getH();
    x=rz*k;
    oscfun = (1.0+tocm->splineBias(zuse))*sin(x)/x;
    oscfun -= 2.0*((3.0/pow(x,3.0)-1.0/x)*sin(x)-3.0*cos(x)/x/x)/3.0;
  }
  oscfun*=c->growthFac(zuse)/c->growthFac(z);

  //get Radio differential flux per atom
  dJdz=getDJradioDz(zuse);

  wk=dJdz*oscfun;
 
  return wk;
  //  return dJdz; //if return dldz then should get wk=1
}

//
double intWindowRadioK(double zuse)
{
  return setWindowRadioK(1.0,zuse,0.0,NULL,NULL,NULL,0);
}


double TwentyOneCM::getBetaTbRadio(double z, double tk, double xi, double lyaflux)
{
  double xa,xc,xtot,xtott;
  double xcHH,xcEH;
  double Tcmb,Tb,Ts,tau;
  double beta;

  Tcmb=TCMB*(1.0+z);
  xa=getXAlpha(z,tk,xi,lyaflux);

 
  xcEH= getXCollEH(z,tk,xi);
  xcHH= getXCollHH(z,tk,xi);
  xc=xcEH+xcHH;

  xtot=xa+xc;
  xtott=xtot*(1.0+xtot);

  Tb=tBrightGen(z,tk,xi,lyaflux);

  Ts=tSpinGen(z,tk,xi,lyaflux);
  tau=tau21CM(z,Ts,xi);

  beta=tau*Ts/(1.0+xtot)/(1.0+z);//This term not zero even if Tb=0

  return beta;
}

double TwentyOneCM::getBetaTbRadioBack(double z, double tk, double xi, double lyaflux)
{
  double Tcmb,Tb,Ts,tau;
  double beta;

  Tcmb=TCMB*(1.0+z);

  Tb=tBrightGen(z,tk,xi,lyaflux);

  Ts=tSpinGen(z,tk,xi,lyaflux);
  tau=tau21CM(z,Ts,xi);

  beta=tau*Tcmb/(1.0+z);//This term not zero even if Tb=0

  return beta;
}

////////////////////////////////////////////////////////////////////
//  Post-reionization 21 cm signal
///////////////////////////////////////////////////////////////////
// Based upon Wyithe and Loeb (2007) 

// Scale dependent factor connecting matter power spectrum to 
// 21 cm power spectrum after reionization
//
// Units: bias   K
double TwentyOneCM::postReionBias21CM(double z, double R)
{
  double bias(1.0e-3);

  return bias;
}

//calculate power spectrum 
//
double TwentyOneCM::postReionPowerSpectrum(double z, double k, double powerDD)
{
  double bias, Dz, ps;
  double step(2.0); //not tracking this exactly, so be crude
  if(z<a->getZReion()+step){
    bias=postReionBias21CM(z,2.0*PI/k);
  }else{
    bias=0.0;
  }
  ps=powerDD*bias*bias;

  return ps;
}
/////////////////////////////////////////////////////////////////
// Calculating full power spectrum
/////////////////////////////////////////////////////////////////

// calculate full power spectrum from components
// angle averaged version
//
double TwentyOneCM::fullPowerSpectrum(double z, double k, double powerFZH, double powerDD, double powerXD, double betaD, double betaX, double betaL, double betaT, double betaV, double gt, double wa, double tk, double xi, double rF, double rT, double powerPostTB)
{
  double powerTB;
  double temp,useT,useL;
  //  double rT, rTcut, rF, rFcut;
  double rTcut, rFcut;
  double xH;
  //  static double rFsave(0.0);
  //static double zsave(-1.0);

  xH=1.0-xi;
  
  powerTB=powerFZH*betaD*betaD;  // FZH contribution
  //slight fudge here with beta coefficients, but betaD~betaX
  
  useT=betaT*gt;
  useL=betaL*wa;
  
  temp=2.0*betaD*useT+useT*useT;
  temp+=2.0*betaD*useL+useL*useL;
  temp+=2.0*useL*useT;
  powerTB+=powerDD*temp;  //DD contribution
  
  temp=2.0*useT+2.0*useL;
  temp*=betaX;
  powerTB+=powerXD*temp;  //XD contribution

  powerTB+=2.0/3.0*betaV*(powerDD*(betaD+useL+useT)+powerXD*betaX);    //angle averaged mu^2 part
  
  powerTB+=powerDD/5.0*betaV*betaV;   // angle averaged mu^4 part

  powerTB*=xH*xH;  //put in missing x_H^2 factor

  //add in post-recombination power spectrum
  powerTB+=powerPostTB;
  
  //apply cutoffs
  //SHOULD ONLY APPLY CUTOFFS TO BARYONIC COMPOENT NOT TO DARK MATTER
  rFcut=pow(1.0+k*k*rF*rF,-2.0); 
  rTcut=exp(-k*k*rT*rT);
  
  powerTB*=rFcut*rTcut;
  return powerTB;
}

//////////////////////////////////////////////////////////
// Transition redshift between lya and collisional coupling
//////////////////////////////////////////////////////////
double TwentyOneCM::transitionZMeanLya()
{
  double ztrans;
  double z1(15.0);
  double z2(40.0);
  double tol(1.0e-4);

  setDummyTransitionZMeanLya(0.0,a,this,1);
  zbracSimp(dummyTransitionZMeanLya,z1,z2);
  ztrans=zriddrSimp(dummyTransitionZMeanLya,z1,z2,tol);

  return ztrans;
}

double setDummyTransitionZMeanLya(double z, Astrophysics *a1, TwentyOneCM *tocm1, int iflag)
{
  static TwentyOneCM *tocm;
  static Astrophysics *a;
  double tk,xe,lya;
  double xa,xc;

  if(iflag==1){
    a=a1;
    tocm=tocm1;
    return 0.0;
  }

  tk=a->getTK(z);
  xe=a->getXE(z);
  lya=a->lyaFlux(z);
  xc=tocm->getXColl(z,tk,xe);
  xa=tocm->getXAlpha(z,tk,xe,lya);

  return xa-xc;
}

double dummyTransitionZMeanLya(double z)
{
  return setDummyTransitionZMeanLya(z,NULL,NULL,0);
}
//////////////////////////////////////////////////////////
// Transition redshift when T_K=T_CMB
//////////////////////////////////////////////////////////
double TwentyOneCM::transitionZMeanT()
{
  double ztrans;
  double z1(5.0);
  double z2(40.0);
  double tol(1.0e-4);

  setDummyTransitionZMeanT(0.0,a,this,1);
  zbracSimp(dummyTransitionZMeanT,z1,z2);
  ztrans=zriddrSimp(dummyTransitionZMeanT,z1,z2,tol);

  return ztrans;
}

double setDummyTransitionZMeanT(double z, Astrophysics *a1, TwentyOneCM *tocm1, int iflag)
{
  static TwentyOneCM *tocm;
  static Astrophysics *a;
  double tk,tcmb;

  if(iflag==1){
    a=a1;
    tocm=tocm1;
    return 0.0;
  }

  tk=a->getTK(z);
  tcmb=TCMB*(1.0+z);

  return tk-tcmb;
}

double dummyTransitionZMeanT(double z)
{
  return setDummyTransitionZMeanT(z,NULL,NULL,0);
}
//////////////////////////////////////////////////////////
// Transition redshift between lya and collisional coupling
//////////////////////////////////////////////////////////
double TwentyOneCM::transitionZLya(double k)
{
  double ztrans;
  double z1(16.0);
  double z2(35.0);
  double tol(1.0e-4);

  Spline wLya;
  Spline wDens;
  double *zN, *flucA, *flucD;
  double step(0.5);
  int i,nz;
  double *betaN;
  double tk,xe,lya;
  double wa;

  z1=transitionZMeanLya();

  nz=(int)((z2-z1)/step);
  nz++;
  step=(z2-z1)/(double)(nz-1);

  zN=dvector(1,nz);
  flucA=dvector(1,nz);
  flucD=dvector(1,nz);
  betaN=dvector(1,5);

  for(i=1;i<=nz;i++){
    zN[i]=z1+(double)(i-1)*step;
    //cout<<zN[i]<<endl;
    tk=a->getTK(zN[i]);
    xe=a->getXE(zN[i]);
    lya=a->lyaFlux(zN[i]);
    wa=getWindowXray(zN[i],k);
    getBeta(zN[i],tk,xe,lya,betaN);
    flucA[i]=wa*betaN[3];
    flucD[i]=betaN[1];
    cout<<zN[i]<<"\t"<<flucA[i]<<"\t"<<flucD[i]<<endl;
  }

  wLya.setSplineSP(nz,zN,flucA);
  wDens.setSplineSP(nz,zN,flucD);

  setDummyTransitionZLya(0.0,&wLya,&wDens,1);
  zbracSimp(dummyTransitionZLya,z1,z2);
  ztrans=zriddrSimp(dummyTransitionZLya,z1,z2,tol);

  return ztrans;
}

double setDummyTransitionZLya(double z, Spline *wLya1, Spline *wDens1, int iflag)
{
  static Spline *wLya;
  static Spline *wDens;

  if(iflag==1){
    wLya=wLya1;
    wDens=wDens1;
    return 0.0;
  }

  return wLya->returnValue(z)-wDens->returnValue(z);
}

double dummyTransitionZLya(double z)
{
  return setDummyTransitionZLya(z,NULL,NULL,0);
}
//////////////////////////////////////////////////////////
// Transition redshift between temp and density fluctuations
//////////////////////////////////////////////////////////
double TwentyOneCM::transitionZT(double k)
{
  double ztrans;
  double z1(10.0);
  double z2(25.0);
  double tol(1.0e-4);

  Spline wTemp;
  Spline wDens;
  double *zN, *flucT, *flucD;
  double step(0.5);
  int i,nz;
  double *betaN;
  double tk,xe,lya;
  double **gtN;

  z1=transitionZMeanT()+0.1;
  //  z2=transitionZMeanLya();

  nz=(int)((z2-z1)/step);
  nz++;
  step=(z2-z1)/(double)(nz-1);

  zN=dvector(1,nz);
  flucT=dvector(1,nz);
  flucD=dvector(1,nz);
  betaN=dvector(1,5);
  gtN=dmatrix(1,4,1,nz);

  //set up redshift values
  for(i=1;i<=nz;i++){
    zN[i]=z1+(double)(i-1)*step;
    cout<<zN[i]<<endl;
  }
  a->getTK(z1);  //make sure global hisoty initialised
  //calculate temperature fluctuations
  getGTN(zN,k,gtN,nz);
  cout<<"done"<<endl;

  //establish splines for temp and density fluctuations
  for(i=1;i<=nz;i++){
    tk=a->getTK(zN[i]);
    xe=a->getXE(zN[i]);
    lya=a->lyaFlux(zN[i]);
    getBeta(zN[i],tk,xe,lya,betaN);
    flucT[i]=gtN[1][i]*betaN[4];
    flucD[i]=betaN[1];
    cout<<zN[i]<<"\t"<<gtN[1][i]<<"\t"<<flucT[i]<<"\t"<<flucD[i]<<endl;
    
  }

  wTemp.setSplineSP(nz,zN,flucT);
  wDens.setSplineSP(nz,zN,flucD);

  setDummyTransitionZT(0.0,&wTemp,&wDens,1);
  
  zbracSimp(dummyTransitionZT,z1,z2);
  ztrans=zriddrSimp(dummyTransitionZT,z1,z2,tol);

  return ztrans;
}

double setDummyTransitionZT(double z, Spline *wTemp1, Spline *wDens1, int iflag)
{
  static Spline *wTemp;
  static Spline *wDens;

  if(iflag==1){
    wTemp=wTemp1;
    wDens=wDens1;
    return 0.0;
  }

  return fabs(wTemp->returnValue(z))-wDens->returnValue(z);
}

double dummyTransitionZT(double z)
{
  return setDummyTransitionZT(z,NULL,NULL,0);
}


///////////////////////////////////////////////////////////
// Post-reionization 21 cm fluctuations
///////////////////////////////////////////////////////////
// Wyithe and Loeb calculations (WL2007)

double TwentyOneCM::tBrightPR(double z, double delta, double R)
{
	double tb, Qi, Fm;
	double deltai;

	tb=22.0e-3*sqrt((1.0+z)/7.5);

	Qi=a->getXI(z);   // not sure this is correct thing to use here
	Fm=fillingMass(deltai,delta,R);
	tb*=1.0-Qi*Fm;
	tb*=1.0+4.0/3.0*delta;
	return tb;
}

//integrate equation (2) of WL2007 to get F_M(\Delta_i)
double TwentyOneCM::fillingMass(double deltai, double delta, double R)
{
	double Fm;


	return Fm;
}

double setFillingMassInt(double z, double deltai1, double delta1, double zobs1, Cosmology *c1, Astrophysics *a1, int iflag)
{
	static Cosmology *c;
	static Astrophysics *a;
	static double delta;
	static double zobs;
	double dFm;

	if(iflag==1){
	a=a1;
	c=c1;
	delta=delta1;
	zobs=zobs1;

	return 0.0;
	}

	return dFm;
}

double fillingMassInt(double z)
{
	return setFillingMassInt(0.0,0.0,0.0,0.0,NULL,NULL,0);
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
//  Calculation of 21cm pdf
///////////////////////////////////////////////////////////

//
// Calculate p(\psi) at specified redshift z, given zeta and pixel size mpix
//
double TwentyOneCM::getPsiPDF(double psi, double z, double zeta, double mpix)
{
   double pdf(0.0);

   

   return pdf;
}

// Use collapse fraction in an overdense region to get ionized fraction
// (3) from FZH04
double TwentyOneCM::XIfromOverdensity(double delta, double z, double zeta, double mpix)
{
   double dsdm;
   double sigmin, sigpix;
   double deltac,deltam;
   double erfc, xarg;
   double xi;

   deltac=c->delCrit0(z)/c->growthFac(z);
   deltam=delta/c->growthFac(z);
   sigmin=sigm(coolMass(c,z),dsdm,c,1);
   sigpix=sigm(mpix,dsdm,c,1);
   //cout<<deltam<<"\t"<<deltac<<"\t"<<sigmin<<"\t"<<sigpix<<endl;
   xarg=(deltac-deltam)/sqrt(2.0*(sigmin*sigmin-sigpix*sigpix));
   erfc=1.0-gsl_sf_erf(xarg);

   xi=zeta*erfc;

   if(xi>1.0) return 1.0;

   return xi;
}

//get psi given the overdensity
//
double TwentyOneCM::psiFromDelta(double delta, double z, double zeta, double mpix)
{
   double xh, psi;

   xh=1.0-XIfromOverdensity(delta,z,zeta,mpix);
   psi=(1.0+delta)*xh;
   return psi;
}

//get psi given the overdensity
//
double TwentyOneCM::derivPsiByDelta(double delta, double z, double zeta, double mpix)
{
   double dpsi;
   double psip, psim;
   double ddelta(0.01);

   ddelta=delta*ddelta;  //fractional may be better, but round off concerns

   psip=psiFromDelta(delta+ddelta,z,zeta,mpix);
   psim=psiFromDelta(delta-ddelta,z,zeta,mpix);
   dpsi=(psip-psim)/(2.0*ddelta);

   return dpsi;
}

//get the two values of delta that give a value of psi
vector<double> TwentyOneCM::deltaFromPsi(double psi, double z, double zeta, double mpix)
{
   double root1(-1.0);
   double root2(5.0);
   double psi1,psi2,psimax;
   double delta;
   double dstep(0.01);
   double psit;
   double delta1, delta2, deltamax;
   double tol(1.0e-3);
   int count;
   vector<double> roots;

   //test limits
   psi1=psiFromDelta(root1,z,zeta,mpix);
   psi2=psiFromDelta(root2,z,zeta,mpix);

   if(psi1>psi || psi2>psi){
      cout<<"limits insufficient"<<endl;
      return roots;
   }

   //get crude estimate of midpoint
   psimax=0.0;
   count=0;
   while(psimax<psi && count<2){
      psit=0.0;
      delta=root1;
      while(delta<root2){
         psit=psiFromDelta(delta,z,zeta,mpix);
         if(psit>psimax){
            psimax=psit;
            deltamax=delta;
         }
         delta+=dstep;
      }

      //cout<<"psimax="<<psimax<<endl;
      dstep/=10.0;
      count++;
   }

   if(psi>psimax){
      cout<<"input value of psi is unreachable"<<endl;
      return roots;
   }

   //can now use psimax for two rootfinding exercises
   setGetPsi(delta,z,zeta,mpix,this,1);
   delta1=zriddrConst(dummyGetPsi,psi,root1,deltamax,tol);
   delta2=zriddrConst(dummyGetPsi,psi,deltamax,root2,tol);
   //cout<<"roots="<<delta1<<"\t"<<delta2<<endl;
   //cout<<"values="<<psiFromDelta(delta1,z,zeta,mpix)<<"\t"<<psiFromDelta(delta2,z,zeta,mpix)<<endl;

   roots.clear();
   roots.push_back(delta1);
   roots.push_back(delta2);
   return roots;
}

double setGetPsi(double delta, double z1, double zeta1, double mpix1, TwentyOneCM *tocm1, int iflag)
{
   double psi;
   static double z;
   static double zeta;
   static double mpix;
   static TwentyOneCM *tocm;
   
   if(iflag){
      z=z1;
      zeta=zeta1;
      mpix=mpix1;
      tocm=tocm1;
      return 0.0;
   }

   psi=tocm->psiFromDelta(delta,z,zeta,mpix);
   return psi;
}

double dummyGetPsi(double delta)
{
   return setGetPsi(delta,0.0,0.0,0.0,NULL,0);
}

double TwentyOneCM::probDelta(double delta, double z, double mpix)
{
   double pdf;
   double sigpix;
   double dsdm;
   double deltam;

   sigpix=sigm(mpix,dsdm,c,1);
   deltam=delta/c->growthFac(z);

   pdf=exp(-deltam*deltam/sigpix/sigpix/2.0)/sqrt(2.0*PI)/sigpix;

   return pdf;
}

//get the raw PDF ignoring the correction for pixels that were ionized 
//externally
double TwentyOneCM::getBrightnessPDF(double psi, double z, double zeta, double mpix)
{
   double pdf(0.0);
   double delta, dPsi, pdelta;
   vector<double> roots;
   double temp;
   roots=deltaFromPsi(psi,z,zeta,mpix);

   for(int i=0;i<roots.size();i++){
      delta=roots[i];
      dPsi=fabs(derivPsiByDelta(delta,z,zeta,mpix));
      pdelta=probDelta(delta,z,mpix);
      temp=pdelta/dPsi;
      cout<<"pdf add="<<delta<<"\t"<<temp<<endl;
      pdf+=temp;
   }
   return pdf/c->growthFac(z);  //unsure on why factor of 1/D(z) here?
}

// Calculate the probability that a pixel was ionized externally
//WARNING: THIS CODE RUNS BUT HAS NOT BEEN VALIDATED - OUTPUT SHOULD NOT BE
//TRUSTED
double TwentyOneCM::getProbPixExternalIon(double siguse, double delta, double z, double zeta, double mpix)
{
   double pion;
   double B,delta0;
   double sigpix, xarg;
   double sigpix2, siguse2;
   double *BB;
   double dsdm;
   BB=dvector(1,3);
   Reionization rz(c,a);

   sigpix=sigm(mpix,dsdm,c,1);

   rz.barrierBub(mpix,z,BB);
   rz.setZeta(zeta);
   B=BB[3];
   // B/=c->growthFac(z); 
   delta0=delta/c->growthFac(z);
   //cout<<"B="<<B<<"\t"<<delta0<<endl;  

   sigpix2=sigpix*sigpix;
   siguse2=siguse*siguse;

   pion=B/sqrt(2.0*PI)*sigpix;
   pion/=siguse*siguse2*sqrt(sigpix2-siguse2);
   xarg=pow(delta0*siguse2-B*sigpix2,2.0)/2.0/siguse2/sigpix2/(sigpix2-siguse2);
   pion*=exp(-xarg);

   free_dvector(BB,1,3);

   //cout<<siguse<<"\t"<<sigpix<<"\t"<<pion<<endl;

   return pion;
}

//integrate the kernel to get the probability that pixel was externally ionized
double TwentyOneCM::evaluateCondProb(double delta, double z, double zeta, double mpix)
{
   double prob;
   double dsdm;
   double lsigmin, lsigmax;
   double tol(1.0e-3);

   //log integral
   lsigmin=-1.0;
   lsigmax=log(0.99999*sigm(mpix,dsdm,c,1));

   setCondProbKernel(0.0,delta,z,zeta,mpix,this,1);

   prob=qromb(getCondProbKernel,lsigmin,lsigmax,tol);

   return prob;
}

double getCondProbKernel(double lsigm)
{
   return setCondProbKernel(lsigm,0.0,0.0,0.0,0.0,NULL,0);
}

double setCondProbKernel(double lsigm, double delta1, double z1, double zeta1, double mpix1, TwentyOneCM *tocm1, int iflag)
{
   double prob, lprob, sigm;
   static double delta;
   static double z;
   static double zeta;
   static double mpix;
   static TwentyOneCM *tocm;

   if(iflag==1){
      delta=delta1;
      z=z1;
      zeta=zeta1;
      mpix=mpix1;
      tocm=tocm1;
      return 0.0;
   }

   sigm=exp(lsigm);
   prob=tocm->getProbPixExternalIon(sigm,delta,z,zeta,mpix);
   //cout<<sigm<<"\t"<<prob<<endl;
   lprob=2.0*sigm*sigm*prob;  // convert d\sigma^2 to 
   return prob;
}

///////////////////////////////////////////////////////////

