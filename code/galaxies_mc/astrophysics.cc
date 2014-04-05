/* astrophysics.cc
 * Contains functions for calculating useful spectroscopic quantites
 * and quantities involved in reionization
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "astrophysics.h"
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

Astrophysics::Astrophysics(Cosmology *c1, int popflag_in, int xrayflag_in, int lyaxray_in, double fxray_in)
{
  // cout <<"Atomic constructor has been called." <<endl;
  c=c1;
  popflag=popflag_in;   //Star type:   0: Pop II 1: Pop III
  xrayflag=xrayflag_in; //xray source: 1: SB  2: SNR  3: mini-quasar
  lyaxray=lyaxray_in;   //lya source:  0: star 1: star+xray 2: xray
  //sourceflag;
  FXRAY=fxray_in;

  if(popflag==1){
    //Pop III stars
    FSTAR=0.01;
    FESC=0.1;
    NION=30000.0;
    NLYA=3030.0;
  }else{
    //Pop II stars
    FSTAR=0.1;
    FESC=0.1;
    NION=4000.0;
    NLYA=6950.0;
  }

  //low tau case
  

  //High tau case
  //   NION=30000.0;
  //FESC=0.1;
  //FSTAR=0.4;

  //cout<<"modified for NB07"<<endl;
  //  FSTAR=0.17;

  if(popflag==1){cout<<"popIII"<<endl;}
  else{cout<<"PopII"<<endl;}

  if(xrayflag==1){cout<<"Starburst"<<endl;}
  else if(xrayflag==2){cout<<"Super nova remnant"<<endl;}
  else if(xrayflag==3){cout<<"Miniquasar"<<endl;}

  if(lyaxray==0){cout<<"stellar lya"<<endl;}
  else if(lyaxray==1){cout<<"stellar+xray lya"<<endl;}
  else if(lyaxray==2){cout<<"xray lya"<<endl;}

  if(MASSFCN==0){cout<<"PS"<<endl;}
  else if(MASSFCN==1){cout<<"ST"<<endl;}
  else if(MASSFCN==2){cout<<"Jenkins"<<endl;}

  ZREION=5.0;
  STARBURSTAGE=1.0e7;
  //cout << "Atomic Constructor called successfully." <<endl; 
  cout<<FSTAR<<"\t"<<FESC<<"\t"<<NION<<"\t"<<NLYA<<"\t"<<FXRAY<<endl;
  globalHistoryFlag=0;
}

Astrophysics::~Astrophysics()
{
  // cout <<"Atomic destructor has been called." <<endl;

}
/////////////////////////////////////////////////////////////////////
// Initialisation routine
/////////////////////////////////////////////////////////////////////
void Astrophysics::initAstrophysics(double fstar, double fesc, double nion, double fx, double nlya, int popflag1, int xrayflag1, int lyaxray1, int sourceflag1, double starburst)
{
  double step(-0.1);
  int cut(-1);

  if(fstar>step) FSTAR=fstar;
  if(fesc>step) FESC=fesc;
  if(nion>step) NION=nion;
  if(fx>step) FXRAY=fx;
  if(nlya>step) NLYA=nlya;
  if(popflag1>cut) popflag=popflag1;
  if(xrayflag1>cut) xrayflag=xrayflag1;
  if(lyaxray1>cut) lyaxray=lyaxray1;
  if(sourceflag>cut) sourceflag=sourceflag1;
  if(starburst>step) STARBURSTAGE=starburst;

  //display choices
  if(popflag==1){cout<<"popIII"<<endl;}
  else{cout<<"PopII"<<endl;}

  if(xrayflag==1){cout<<"Starburst"<<endl;}
  else if(xrayflag==2){cout<<"Super nova remnant"<<endl;}
  else if(xrayflag==3){cout<<"Miniquasar"<<endl;}

  if(lyaxray==0){cout<<"stellar lya"<<endl;}
  else if(lyaxray==1){cout<<"stellar+xray lya"<<endl;}
  else if(lyaxray==2){cout<<"xray lya"<<endl;}

  if(MASSFCN==0){cout<<"PS"<<endl;}
  else if(MASSFCN==1){cout<<"ST"<<endl;}
  else if(MASSFCN==2){cout<<"Jenkins"<<endl;}

  cout<<FSTAR<<"\t"<<FESC<<"\t"<<NION<<"\t"<<NLYA<<"\t"<<FXRAY<<endl;

  globalHistoryFlag=0;
}

void Astrophysics::setFstar(double fstar)
{
  initAstrophysics(fstar,-1.0,-1.0,-1.0,-1.0,-1,-1,-1,-1,-1.0);
}

void Astrophysics::setFX(double fx)
{
  initAstrophysics(-1.0,-1.0,-1.0,fx,-1.0,-1,-1,-1,-1,-1.0);
}

void Astrophysics::setNION(double nion)
{
  initAstrophysics(-1.0,-1.0,nion,-1.0,-1.0,-1,-1,-1,-1,-1.0);
}


void Astrophysics::setNLYA(double nlya)
{
  initAstrophysics(-1.0,-1.0,-1.0,-1.0,nlya,-1,-1,-1,-1,-1.0);
}

void Astrophysics::setCosmology(Cosmology *c1)
{
  c=c1;
}

/////////////////////////////////////////////////////////////////////
// Member functions
/////////////////////////////////////////////////////////////////////

double Astrophysics::getFSTAR(void)
{
  return FSTAR;
}

double Astrophysics::getFESC(void)
{
  return FESC;
}

double Astrophysics::getNION(void)
{
  return NION;
}

double Astrophysics::getNLYA(void)
{
  return NLYA;
}

int Astrophysics::getPopflag(void)
{
  return popflag;
}

double Astrophysics::getZeta(void)
{
  return FSTAR*FESC*NION;
}

int Astrophysics::getLyaXray(void)
{
  return lyaxray;
}

double Astrophysics::getFXRAY(void)
{
  return FXRAY;
}

double Astrophysics::clumpingIntMEHR(double ldeltai)
{
  return exp(MEHRclumping.returnValue(ldeltai));
}

double Astrophysics::volumeIntMEHR(double ldeltai)
{
  return exp(MEHRvolume.returnValue(ldeltai));
}

double Astrophysics::getZReion()
{
  return ZREION;
}


/*********************************************************************
 **************** Utility Functions      *****************************
 ********************************************************************/
// Calculates the transition probability A_{ul} between an upper 
// state of hydrogen  u=(n,l,j) and lower state l=(n',l',j').  
// 
double Astrophysics::transitionA(double n, double l, double j, double np, double lp,
			    double jp)
{
  double tranA(1.0);
  double lmax;
  double coupling;

  coupling=gsl_sf_coupling_6j(int(2.0*l),int(2.0*j),int(2.0*0.5),int(2.0*jp),int(2.0*lp),int(2.0*1.0));

  if(coupling == 0.0){
    return 0.0;
  }

  lmax=fmax(fabs(l),fabs(lp));
  tranA = lmax*(2.0*jp+1.0)*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS;
  tranA *= BOHRA0*BOHRA0*coupling*coupling;
  tranA *= pow(hydrogenMatrix(n,l,np,lp),2.0);
  tranA *= 64.0*PI*PI*PI*PI/3.0/PLANCKCONSTANT;
  tranA /= pow(hydrogenWave(n,np),3.0);
  
  return tranA;
}

//Need to modify below so that it works for (n,l) to (np,lp) transitions.

// Calculates the absorption oscillator strength of the transition  between
// an upper state of hydrogen  u=(n,l,j) and lower state l=(n',l',j')
double Astrophysics::transitionF(double n, double l, double j, double np, double lp,
			    double jp)
{
  double transF;

  transF = ELECTRONMASS*pow(SPEEDOFLIGHT_CGS,3.0);
  transF/= 8.0*PI*PI*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS*pow(hydrogenFreq(n,np),2.0);
  transF *= transitionA(n,l,j,np,lp,jp);
  transF *= (2.0*j+1.0)/(2.0*jp+1.0);

  return transF;
} 

// Calculates the frequency of the transition between the upper level n
// and lower level np in the hydrogen atom.
double Astrophysics::hydrogenFreq(double n, double np)
{
  double hf(1.0);
  
  hf=RYDBERG_CGS*(1.0/np/np - 1.0/n/n)/PLANCKCONSTANT;

  return hf;
}

// Calculates the wavelength of the transition between the upper level n
// and lower level np in the hydrogen atom.
double Astrophysics::hydrogenWave(double n, double np)
{
  double hl(SPEEDOFLIGHT_CGS);
  
  hl /= RYDBERG_CGS*(1.0/np/np - 1.0/n/n)/PLANCKCONSTANT;
  
  return hl;
}

// Calculates the Energy of the transition between the upper level n
// and lower level np in the hydrogen atom.
double Astrophysics::hydrogenEnergy(double n, double np)
{
  double he(RYDBERG_CGS);
  
  he *= (1.0/np/np - 1.0/n/n);
  
  return he;
}


// Calculates the transition strength S_{ul} between an upper 
// state of hydrogen  u=(n,l,j) and lower state l=(n',l',j').  
//By convention the 
// 
double Astrophysics::transitionS(double n, double l, double j, double np, double lp,
			    double jp)
{
  double tranS(1.0);
  double lmax;
  double coupling;

  coupling=gsl_sf_coupling_6j(int(2.0*l),int(2.0*j),int(2.0*0.5),int(2.0*jp),int(2.0*lp),int(2.0*1.0));

  if(coupling == 0.0){
    return 0.0;
  }

  lmax=fmax(fabs(l),fabs(lp));
  tranS = lmax*(2.0*j+1.0)*(2.0*jp+1.0)*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS;
  tranS *= BOHRA0*BOHRA0*coupling*coupling;
  tranS *=pow(hydrogenMatrix(n,l,np,lp),2.0);
  
  return tranS;
}

// Calculates the total decay rate from a state u=(n,l,j).
//
// Loop structure could probably be made more elegant
double Astrophysics::transitionATotal(double n, double l, double j)
{
  double Atotal(0.0);
  double m,jp;
  
  // For l>0 case handle transitions which have Delta l = -1
  if(l > 0.0){
    m=l;
    while(m<=n-1.0){
      jp=fabs(l-1.5);
      while(jp<=l-0.5){
	Atotal+=transitionA(n,l,j,m,l-1.0,jp);
	jp += 1.0;
      }
      	m += 1.0;
    }
  }
  // For all cases handle transitions with Delta l = +1
  m=l+2.0;
  while(m<=n-1.0){
    jp=fabs(l+0.5);
    while(jp<=l+1.5){
      Atotal+=transitionA(n,l,j,m,l+1.0,jp);
      jp += 1.0;
    }
    m += 1.0;
  } 
 
  return Atotal;

}

// Calculates the probability P_ul that a photon in the state u=(n,l,j) will
// make a transition to the state l=(n',l',j')
  double Astrophysics::transitionP(double n, double l, double j, double np,
			      double lp, double jp)
{
  double transP(1.0);

  transP=transitionA(n,l,j,np,lp,jp)/transitionATotal(n,l,j);

  return transP;
}


// Calculates the matrix element R^{n',l'}_{n,l}/a_0 for the hydrogen 
// atom.  Here (n',l') describe the lower level and (n,l) the upper level.
// a_0 is the Bohr radius.
// Reference: Rudnick (1935) Physical Review 48, 807
double Astrophysics::hydrogenMatrix(double n, double l, double np, double lp)
{
  double matR;

  matR=factorial(n+l)/factorial(np+lp)/factorial(n-l-1.0)/factorial(np-lp-1.0);
  matR=sqrt(matR);
  matR *= pow(n-np,n-np-1.0)/pow(n+np,n+np+1.0);
  matR *= pow(2.0,l+lp+4.0)*pow(n,lp+3.0)*pow(np,l+3.0);
  matR *= polyR(n,l,np,lp);
  return matR;
}

// Calculates the polynomial required by hydrogenMatrix 
double Astrophysics::polyR(double n, double l, double np, double lp)
{
  double pR;
  double lmax,v,nr,r;

  lmax=fmax(fabs(l),fabs(lp));
  r=np-lmax-1.0;
  nr=n-l-1.0;
  v=-4.0*n*np/pow(n-np,2.0);

  if(fabs(lp-l+1.0)<1.0e-5){
    //This is the delta l = -1 case
    pR = (n+np)*term2F1(-r,-n+lmax+1.0,2*lmax+1.0,v);
    pR -= (n-np)*term2F1(-r,-n+lmax,2*lmax+1.0,v);
    pR *=pow(-1.0,r)*pow(n-np,2.0*r)*factorial(2.0*lmax+r)/factorial(2.0*lmax);
    pR /= 2.0*np;
    return pR;
  }

   if(fabs(lp-l-1.0)<1.0e-5){
    //This is the delta l = +1 case
    pR = (n+np)*(n-lmax)*term2F1(-r,-n+lmax+1.0,2*lmax+1.0,v);
    pR -= (n-np)*(n+lmax)*term2F1(-r,-n+lmax,2*lmax+1.0,v);
    pR *=pow(-1.0,r)*pow(n-np,2.0*r)*factorial(2.0*lmax+r)/factorial(2.0*lmax);
    pR /= 2.0*n;
    return pR;
  } 

   cout << "Delta l not equal to +- 1 in PolyR!\n" ;
   return 0.0;
}

/*********************************************************************
 ********* Math functions for Hypergeometric Functions   *************
 ********************************************************************/
//Calculates the value of the Pochhammer symbol (a)_k
double Astrophysics::intPochhammer(double a, double k)
{
  double temp(1.0);
  double j(0.0);

  if (k <0.0) {
    cout << "k negative in intPochhammer!\n";
    return 0.0;
  }
  if(k ==0.0) {
    return 1.0;
  }

  while(j <= k-1.0){
    temp *= (a+j);
    j +=1.0;
      }

  return temp;
}

//Calculate the factorial n!
double Astrophysics::factorial(double n)
{
  double temp(1.0);
  double j(1.0);

  if(n<0.0){
    cout << "negative argument to factorial!\n";
    return 0.0;
  }

  while(j<n){
    j += 1.0;
    temp *= j;
  }
  return temp;
}
	

// Calculates the hypergeometric function 2F1(a,b;c,x)
// Requires that a or b be a negative integer so that the series
// terminates
double Astrophysics::term2F1(double a, double b, double c, double x)
{
  double temp(1.0);
  double k(1.0);
  double kmax(0.0);

  // Check to make sure that at least one of a or b are negative
  if( (a>0.0) && (b>0.0)){
    cout << "a and b both positive in term2F1. Series must terminate\n";
    return 0.0;
  }
  // Work out when series terminates
  if( (a>0.0) && (b<0.0)) {
    kmax=-b;
  }
   if( (a<0.0) && (b>0.0)) {
    kmax=-a;
  } 
   if( (a<0.0) && (b<0.0)) {
    kmax=fmin(fabs(a),fabs(b));
  } 

   while( k<= kmax){
     temp += intPochhammer(a,k)*intPochhammer(b,k)*pow(x,k)/intPochhammer(c,k)/factorial(k);
     k +=1.0;
   }
   return temp;
}

//////////////////////////////////////////////////////////////////////////
//  Photo ionization cross-section
//////////////////////////////////////////////////////////////////////////
//Ref: Spitzer, Chap. 5


//Calculate the photoionisation cross-section from the ground state
// for an atom of charge Z as a function of frequency in Hertz
//Units:   nu    Second^{-1}
//         sigma Centimeter^2
double Astrophysics::xsectionPhotoIonise(double nu, double Z)
{
  double sigma(7.91e-18);
  double x,z;
  double nu1(Z*Z*NULYMANLIMIT);
  
  //handle frequencies at or just above edge to avoid divide by zero
  if(nu-nu1<1.0e-4) return 0.0;

  x=nu1/nu;
  z=sqrt(nu1/(nu-nu1));

  sigma*=x*x*x/Z/Z;
  sigma*=gaunt1F(x,z);
  return sigma;
}

//Calculate the gaunt factor for transitions from the n=1 level to the
//continuum for a hydrogen-like atom.
// x=nu1/nu, z=nu1/(nu=nu1)
//Spitzer: p105.  Units: x,z  1
double Astrophysics::gaunt1F(double x, double z)
{
  double g1F;
  double temp;
  
  g1F=8.0*PI*sqrt(3.0)*x;
  temp=exp(-4.0*z*atan2(1.0,z));   //arccot(z)=arctan(1/z)
  temp/=(1.0-exp(-2.0*PI*z));
  g1F*=temp;
  return g1F;
}

//calculate the mean free path for photo ionisation for a photon
//Units:  energy eV 
//        mfp    Mpc
double Astrophysics::mfpPhotoIonise(double energy, double z, double xfree)
{
  double mfp;
  double nu;

  nu=energy/RYDBERG*NULYMANLIMIT;

  mfp=1.0/xsectionPhotoIonise(nu,1.0);
  mfp/=c->nh(z);
  mfp/=3.0856e24;  //convert cm into Mpc
  return mfp;
}

//General fits to photoionization cross-sections
//taken from Verner et al. (1996)
//
//Units : xsection  cm^2
double Astrophysics::xsectionPhotoIoniseGen(double nu, int z, int n)
{
  double xsection,E;
  double Eth(0.0), Emax(0.0), E0(0.0);
  double sigma0(0.0),ya(0.0),yw(0.0),y0(0.0),y1(0.0),P(0.0);

  double x,y,F;
  
  //Find fitting values
  if(z==1 && n==1){  //HI
    Eth=1.36e1;
    Emax=5.0e4;
    E0=4.298e-1;
    sigma0=5.475e4;
    ya=3.288e1;
    P=2.963e0;
    yw=0.0;
    y0=0.0;
    y1=0.0;
  }else if(z==2 && n==1){   //HeI
    Eth=2.459e1;
    Emax=5.0e4;
    E0=1.361e1;
    sigma0=9.492e2;
    ya=1.469e0;
    P=3.188e0;
    yw=2.039e0;
    y0=4.434e-1;
    y1=2.136e0;
  }else if(z==2 && n==2){   //HeII
    Eth=5.442e1;
    Emax=5.0e4;
    E0=1.720e0;
    sigma0=1.369e4;
    ya=3.288e1;
    P=2.963e0;
    yw=0.0;
    y0=0.0;
    y1=0.0;
  }else{
    cout<<"no data for z>2 in xsectionPI"<<endl;
    return 0.0;
  }

  E=nu/NULYMANLIMIT*RYDBERG;
  if(E<Eth) return 0.0;  // E must be greater than threshold
  x=E/E0-y0;
  y=sqrt(x*x+y1*y1);
  F=(x-1.0)*(x-1.0)+yw*yw;
  F*=pow(y,0.5*P-5.5);
  F*=pow(1.0+sqrt(y/ya),-P);

  xsection=sigma0*F*1.0e-18;

  return xsection;
}
    
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//  Black body functions
///////////////////////////////////////////////////////////////////////////
//UNCHECKED

//BB energy density
//Units:  nu Second^-1
//        T  Kelvin
//        u  erg cm^-3 
double Astrophysics::uBlackBody(double nu, double T)
{
  double x,u;

  x=PLANCKCONSTANT*nu/BOLTZK/T;
  u=2.0*PLANCKCONSTANT/SPEEDOFLIGHT_CGS/SPEEDOFLIGHT_CGS/SPEEDOFLIGHT_CGS;
  u*=nu*nu*nu;
  u/=exp(x)-1.0;

  return u;
}


//BB energy flux
//Units:  I_nu erg cm^-2 s^-1 ster^-1
double Astrophysics::fluxBlackBody(double nu, double T)
{
  return SPEEDOFLIGHT_CGS*uBlackBody(nu,T);
}

//BB number flux
//Units: N_nu cm^-2 s^-1 ster^-1
double Astrophysics::nfluxBlackBody(double nu, double T)
{
  double N;
  N=fluxBlackBody(nu,T)/PLANCKCONSTANT/nu;
  return N;
}

//BB number flux
//Units: N_nu cm^-3
double Astrophysics::ndensityBlackBody(double nu, double T)
{
  double N;
  N=nfluxBlackBody(nu,T)*4.0*PI/SPEEDOFLIGHT_CGS;
  return N;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//  Starburst galaxy functions
///////////////////////////////////////////////////////////////////////////
//Total luminosity of a star burst galaxy - in x-ray band
//I take a power law model based on Glover and Brand (2003)
//calibrated for 2-10keV range, but applicable outside that
//Units:  nu  s^-1
//        Ltot   ergs s^-1
//        SFR in Msun yr^-1
double Astrophysics::lumTotStarBurst(double SFR)
{
  double L0(3.4e40);  //normalisation
  double Ltot;

  Ltot=L0*SFR;

  return FXRAY*Ltot;
}


//Specific luminosity of a star burst galaxy - in x-ray band
//I take a power law model based on Glover and Brand (2003)
//calibrated for 2-10keV range, but applicable outside that
//Units:  nu  s^-1
//        Lnu   ergs s^-1 Hz^{-1}
//        SFR in Msun yr^-1
double Astrophysics::lumNuStarBurst(double nu, double SFR)
{
  double L0(3.4e22);  //normalisation in erg s^-1 Hz^-1
  double alpha(1.5);  //spectral index
  double E0(1.0e3);   //pivot energy in eV
  double E,Lnu;
  E=nu/NULYMANLIMIT*RYDBERG;

  Lnu=L0*pow(E0/E,alpha);
  Lnu*=SFR;

  return FXRAY*Lnu;
}

//flux from star burst galaxy as a function of distance
//
//  Units: nu S^-1
//         r  Mpc
//         Jnu  erg cm^-2 s^-1 Hz^-1 ster^-1
double Astrophysics::jNuStarBurst(double nu, double r, double z)
{
  double Jnu;
  double tau;
  double rcm;
  rcm=r*MPC;
  
  tau=rcm*c->nh(z)*xsectionPhotoIonise(nu,1.0);

  Jnu=lumNuStarBurst(nu,1.0)/4.0/PI/rcm/rcm;
  Jnu*=exp(-tau);

  cout<<nu<<"\t"<<r<<"\t"<<tau<<"\t"<<z<<endl;

  return Jnu;
}

//Calculate the specific heating rate for photoionisation of Hydrogen
//from a starburst galaxy
//Units : lambda  erg s^{-1}Hz^{-1}
//         r      Mpc
double Astrophysics::heatingNuStarBurst(double nu, double r, double z)
{
  double lambda;
  double fheat(1.0);

  lambda=4.0*PI*c->nh(z)*xsectionPhotoIonise(nu,1.0);
  lambda*=jNuStarBurst(nu,r,z);
  lambda*=(nu-NULYMANLIMIT)/nu;
  lambda*=fheat;

  return lambda;
}

//Calculate the total heating rate from photoionisation of Hydrogen
//from a starburst galaxy
//Units:  lambda erg s^{-1}Hz^{-1}
//        r  Mpc
double Astrophysics::heatingStarBurst(double r, double z)
{
  double lambda;
  double lEmin,lEmax;
  double tol(1.0e-4);

  lEmin=log(3.0e2);   //0.3-30 keV band
  lEmax=log(3.0e4);

  setLambdaNuSB(0.0,r,z,this,1);
  
  lambda=qromb(getLambdaNuSB,lEmin,lEmax,tol);

  lambda*=RYDBERG_CGS/RYDBERG;      //convert from electron volts into Ergs

  return lambda/PLANCKCONSTANT;
}

double getLambdaNuSB(double lE)
{
  return setLambdaNuSB(lE,0.0,0.0,NULL,0);
}

double setLambdaNuSB(double lE, double r1, double z1, Astrophysics *a1,int iflag)
{
  static double r;
  static double z;
  static Astrophysics *a;
  double E,nu;

  if(iflag==1){
    r=r1;
    z=z1;
    a=a1;
    return 0.0;
  }
  E=exp(lE);
  nu=E/RYDBERG*NULYMANLIMIT;

  return a->heatingNuStarBurst(nu,r,z)*E;
}

///////////////////////////////////////////////////////////////////////////
//  Super nova remnants
///////////////////////////////////////////////////////////////////////////
//Total luminosity of a SNr - in x-ray band
//I take a power law model based on Glover and Brand (2003)
//calibrated to match starburst luminosity in 0.2-30keV range
//Units:  nu  s^-1
//        Ltot   ergs s^-1
//        SFR in Msun yr^-1
double Astrophysics::lumTotSNr(double SFR)
{
  double L0(3.4e40);  //normalisation
  double Ltot;

  Ltot=L0*SFR;

  return FXRAY*Ltot;
}


//Specific luminosity of a SNr - in x-ray band
//I take a power law model based on Glover and Brand (2003)
//calibrated to match total starburst luminosity in 0.2-30 keV range
//Units:  nu  s^-1
//        Lnu   ergs s^-1 Hz^{-1}
//        SFR in Msun yr^-1
double Astrophysics::lumNuSNr(double nu, double SFR)
{
  double L0(2.8e22);  //normalisation in erg s^-1 Hz^-1
  double alpha(1.0);  //spectral index
  double E0(1.0e3);   //pivot energy in eV
  double E,Lnu;
  E=nu/NULYMANLIMIT*RYDBERG;

  Lnu=L0*pow(E0/E,alpha);
  Lnu*=SFR;

  return FXRAY*Lnu;
}

///////////////////////////////////////////////////////////////////////////
//  Mini-Quasar
///////////////////////////////////////////////////////////////////////////
//Total luminosity of a Mini-Quasar - in x-ray band
//I take a power law model based on Glover and Brand (2003)
//calibrated to match starburst luminosity in 0.2-30keV range
//Units:  nu  s^-1
//        Ltot   ergs s^-1
//        SFR in Msun yr^-1
double Astrophysics::lumTotMiniQuasar(double SFR)
{
  double L0(3.4e40);  //normalisation
  double Ltot;

  Ltot=L0*SFR;

  return FXRAY*Ltot;
}


//Specific luminosity of a MiniQuasar - in x-ray band
//I take a power law model based on Glover and Brand (2003)
//calibrated to match total starburst luminosity in 0.2-30 keV range
//Units:  nu  s^-1
//        Lnu   ergs s^-1 Hz^{-1}
//        SFR in Msun yr^-1
double Astrophysics::lumNuMiniQuasar(double nu, double SFR)
{
  double L0(1.4e22);  //normalisation in erg s^-1 Hz^-1
  double alpha(0.5);  //spectral index
  double E0(1.0e3);   //pivot energy in eV
  double E,Lnu;
  E=nu/NULYMANLIMIT*RYDBERG;

  Lnu=L0*pow(E0/E,alpha);
  Lnu*=SFR;

  return FXRAY*Lnu;
}
////////////////////////////////////////////////////////////////////////
//  Calculate Lyman alpha flux at a given redshift
///////////////////////////////////////////////////////////////////////
//calculate mean Lya flux at redshift z
// Reference: Barkana and Loeb (2005) detecting
//Units: lyaflux cm^3 s^-1 Hz^-1 ster^-1
double Astrophysics::lyaFlux(double z)
{
  double lyaflux(0.0);
  double tol(1.0e-4);
  double n(2.0);
  double nmax((double)(lynmax));
  double zmin,zmax;
  setDJalphaDz(z,z,c,this,1);

  zmin=z+getZHII(z);
  n=nmax;
  while(n>=2.0){
    zmax=lymanZMax(z,n);
    lyaflux+=qromb(getDJalphaDz,zmin,zmax,tol);
    zmin=zmax;
    n-=1.0;
  }  

  return lyaflux;
}

// Function to calculate the sourceEmission, eta
// Units: nu s^-1
//        sE cm^-3 s^-1 Hz^-1
// might be better to make this per unit SFR
double Astrophysics::sourceEmission(double nu, double z)
{
  double emission0(1.958e-16);
  double spec_alpha1(1.29);
  double spec_alpha2(0.554);
  double zz;
  double nulymanbeta(hydrogenFreq(3.0,1.0));
  double sE;
  double nlyanorm;
  zz=1.0+z;
  
  if(popflag==0){
    // Pop II spectrum
    emission0=265.4;
    spec_alpha1=0.14;
    spec_alpha2=-6.6589;
    nlyanorm=6950.0;
  }else if(popflag==1){
    //Pop III spectrum
    emission0=1.958e-16;
    spec_alpha1=1.29;
    spec_alpha2=0.554;
    nlyanorm=3030.0;
  }else if(popflag==2){
    //Pop III spectrum - no break
    emission0=1.958e-16;
    spec_alpha1=1.29;
    spec_alpha2=1.29;
    nlyanorm=3030.0;  //NEEDS UPDATING
  }

  sE=emission0;
  sE*=NLYA/nlyanorm;  //rescale emission to get NLYA lya photons per baryon

  if(nu<nulymanbeta){
    sE *= pow(nu,spec_alpha1-1.0);
  }else if(nu>NULYMANLIMIT){
    cout<<"nu bigger than lyman limit in sourceEmission"<<endl;
    sE=0.0;
  }else{
    sE*=pow(nulymanbeta,spec_alpha1-1.0)*pow(nu/nulymanbeta,spec_alpha2-1.0);
  }

  sE *= globalSFR(z)/zz/zz/zz/PROTONMASS;  //comoving SFR rate
  sE/=1.0e9*KPC*KPC*KPC*YEAR/SOLARMASS;

   return sE;
}

// Function to calculate the differential lyn flux on a gas element at 
// redshift z, from a source at zp
double setDJalphaDzDn(double n1, double z1, double zp, Cosmology *c1, Astrophysics *a1,int iflag)
{
  static double n;
  static double z;
  static Cosmology *c;
  static Astrophysics *a;
  double dJdz;
  double nup;

  if(iflag==1){
    n=n1;
    z=z1;
    a=a1;
    c=c1;
    if(zp<3.0) zp=z;  //handle outside of spline values
    //    return 0.0;
  }
  if(zp>a->lymanZMax(z,n)) return 0.0;

  nup=a->hydrogenFreq(n,1.0)*(1.0+zp)/(1.0+z);

  dJdz=SPEEDOFLIGHT_CGS/4.0/PI;
  dJdz *= (1.0+z)*(1.0+z)*a->sourceEmission(nup,zp);
  dJdz /= c->hubbleZ(zp)*UNH*c->getH();
  
  return dJdz;
}

double getDJalphaDzDn(double zp)
{
  return setDJalphaDzDn(0.0,0.0,zp,NULL,NULL,0);
}

//function to calculate the full lya flux at redshift z from redshift zp
double setDJalphaDz(double z1,double zp,Cosmology *c1, Astrophysics *a1, int iflag)
{
  static Cosmology *c;
  static Astrophysics *a;
  static double z;
  double lyaflux(0.0);
  double n(2.0);
  double nmax((double)(lynmax));
  static double xe;

  if(iflag==1){
    c=c1;
    a=a1;
    z=z1;
    xe=a->getXE(z);
    return 0.0;
  }

  //sum over different levels
  if(a->getLyaXray()<2){
    while(n<=nmax){
      setDJalphaDzDn(n,z,0.0,c,a,1);
      lyaflux+=getDJalphaDzDn(zp)*a->lyaRecycle(n);
      n+=1.0;
    }
  }

  if(a->getLyaXray()>0){
    setDJxrayDzDE(10.0,z,0.0,c,a,1);
    getDJalphaXrayDzDE(xe,z,a,c,1);
    lyaflux+=getDJalphaXrayDz(zp);
  }
  return lyaflux;
}

double getDJalphaDz(double zp)
{
  return setDJalphaDz(0.0,zp,NULL,NULL,0);
}


double Astrophysics::lymanZMax(double z, double n)
{
  double zlmax;

  zlmax=(1.0+z)*(1.0-pow(n+1.0,-2.0))/(1.0-pow(n,-2.0));
  zlmax -= 1.0;

  return zlmax;
}

//Return the recycling fraction for photons redshifting into a level n
double Astrophysics::lyaRecycle(double n)
{
  if(blflag==1) return 1.0;

  if(fabs(n-2.0)<1.0e-6) return 1.0;
  if(fabs(n-3.0)<1.0e-6) return 0.0;
  if(fabs(n-4.0)<1.0e-6) return 0.2609;
  if(fabs(n-5.0)<1.0e-6) return 0.3078;
  if(fabs(n-6.0)<1.0e-6) return 0.3259;
  if(fabs(n-7.0)<1.0e-6) return 0.3353;
  if(fabs(n-8.0)<1.0e-6) return 0.3410;
  if(fabs(n-9.0)<1.0e-6) return 0.3448;
  if(fabs(n-10.0)<1.0e-6) return 0.3476;
  if(fabs(n-11.0)<1.0e-6) return 0.3496;
  if(fabs(n-12.0)<1.0e-6) return 0.3512;
  if(fabs(n-13.0)<1.0e-6) return 0.3524;
  if(fabs(n-14.0)<1.0e-6) return 0.3535;
  if(fabs(n-15.0)<1.0e-6) return 0.3543;
  if(fabs(n-16.0)<1.0e-6) return 0.3550;
  if(fabs(n-17.0)<1.0e-6) return 0.3556;
  if(fabs(n-18.0)<1.0e-6) return 0.3561;
  if(fabs(n-19.0)<1.0e-6) return 0.3565;
  if(fabs(n-20.0)<1.0e-6) return 0.3569;
  if(fabs(n-21.0)<1.0e-6) return 0.3572;
  if(fabs(n-22.0)<1.0e-6) return 0.3575;
  if(fabs(n-23.0)<1.0e-6) return 0.3578;
  if(fabs(n-24.0)<1.0e-6) return 0.3580;
  if(fabs(n-25.0)<1.0e-6) return 0.3582;
  if(fabs(n-26.0)<1.0e-6) return 0.3584;
  if(fabs(n-27.0)<1.0e-6) return 0.3586;
  if(fabs(n-28.0)<1.0e-6) return 0.3587;
  if(fabs(n-29.0)<1.0e-6) return 0.3589;
  if(fabs(n-30.0)<1.0e-6) return 0.3590;
  if(n>30.0) return 0.3590;

  cout<<"error in lyaRecycle"<<endl;
  return 0.0;
}

//calculate the delta redshift corresponding to the size of an HII region 
// for a gas element at redshift z
double Astrophysics::getZHII(double z)
{
  double rHII,zHII,dz;
  
  rHII=0.01/pow(z/21.0,2.0);  // comoving size of HII region 10kpc @ z=20
  rHII*=10.0;
  zHII=c->zOfR(rHII,z);
  dz=zHII-z;

  //  return 0.0;
  return dz;
}

////////////////////////////////////////////////////////////////////////
//  Calculate Xray flux at a given redshift
///////////////////////////////////////////////////////////////////////
//calculate mean xray flux at redshift z
//Units: xrayflux cm^3 s^-1 Hz^-1 ster^-1
double Astrophysics::xrayFlux(double z)
{
  double xrayflux(0.0);
  double tol(1.0e-4);
  double zstart;
  double zmax(ZSFR);
  //  double *result;
  //double tk,xi,xe;

  //result=dvector(1,3);
  //getTIGM(z,result);
  //tk=result[1];
  //xi=result[2];   
  //xe=result[2];
  //free_dvector(result,1,3);

  //sum over different levels
  setDJxrayDzDE(20.0,z,0.0,c,this,1);
  zstart=z;
  xrayflux=qromb(getDJxrayDz,zstart,zmax,tol);
  
  return xrayflux;
}

// Function to calculate the sourceEmission, eta for xrays
// Units: E eV
//        sE cm^-3 s^-1 Hz^-1 (comoving)
// might be better to make this per unit SFR
double Astrophysics::sourceEmissionXray(double E, double z)
{
  double zz;
  double sE;
  double SFR,nu;
  zz=1.0+z;

  if(E<XrayEmin || E>XrayEmax) return 0.0;  //limits to xray emission band 

  nu=E/RYDBERG*NULYMANLIMIT;

  SFR=globalSFR(z)/zz/zz/zz;  //comoving SFR rate

  //Choose x-ray source
  if(xrayflag==1){ sE=lumNuStarBurst(nu,SFR);}  
  else if(xrayflag==2){ sE=lumNuSNr(nu,SFR); }
  else if(xrayflag==3){ sE=lumNuMiniQuasar(nu,SFR); }
  else {sE=0.0;}


  sE/=MPC*MPC*MPC; //convert to erg cm^-3 s^-1 Hz-1
  sE/=E/RYDBERG*RYDBERG_CGS;  //convert to cm^-3 s^-1 Hz^-1

  return FXRAY*sE;
}

// Function to calculate the differential xray flux on a gas element at 
// redshift z, from a source at zp
double setDJxrayDzDE(double E, double z1, double zp, Cosmology *c1, Astrophysics *a1,int iflag)
{
  static double z;
  static Cosmology *c;
  static Astrophysics *a;
  double dJdz;
  double Ep;
  double tau;

  if(iflag==1){
    z=z1;
    a=a1;
    c=c1;
    a=a1;
    if(zp<3.0) zp=z;  //handle outside of spline values
    //    return 0.0;
  }

  Ep=E*(1.0+zp)/(1.0+z);

  dJdz=SPEEDOFLIGHT_CGS/4.0/PI;
  dJdz *= (1.0+z)*(1.0+z)*a->sourceEmissionXray(Ep,zp);
  dJdz /= c->hubbleZ(zp)*UNH*c->getH();
  dJdz/=PLANCKCONSTANT;  //working in energy not nu

  //Account for optical depth
  //THIS REALLY SLOWS DOWN MY CODE!!!!
  tau=a->XrayTau(z,zp,Ep);
  dJdz*=exp(-tau);  

  return dJdz;
}


//returns differential xray number flux as function of energy E and
//emission redshift zp
double getDJxrayDzDE(double E, double zp)
{
  return setDJxrayDzDE(E,0.0,zp,NULL,NULL,0);
}

//global variable for double integral in xray calculations
double zpsave;

//return dJdz having integrated over energies
double getDJxrayDz(double zp)
{
  double Emin(RYDBERG);
  double Emax(XrayEmax);
  double dJdz;
  double lEmin,lEmax;
  double tol(1.0e-4);

  lEmin=log(Emin);
  lEmax=log(Emax);

  zpsave=zp;

  dJdz=qromb(getDJxrayDlE,lEmin,lEmax,tol);

  return dJdz;
}

double getDJxrayDlE(double lE)
{
  double E,zp;
  E=exp(lE);
  zp=zpsave;

  return getDJxrayDzDE(E,zp)*E*RYDBERG_CGS/RYDBERG;
}

/////////////////////////////////////////////////////
// X-ray Optical depths
/////////////////////////////////////////////////////////


// Calculate optical depth between redshift z and zp
// for an X-ray of energy E
// ONLY VALID ON SMALL SCALES, IF AT ALL!!!
double Astrophysics::XrayTauSimp(double z, double zp, double E)
{
  double tau,rcm;
  double nu;

  nu=E/RYDBERG*NULYMANLIMIT;
  rcm=c->relativeConformalR(z,zp)*MPC;
  tau=rcm*c->nh(z)*xsectionPhotoIonise(nu,1.0);

  return tau;
}

// Calculate optical depth between redshift z and zp
// for an X-ray emitted at ze with energy Ee
double Astrophysics::XrayTau(double z, double ze, double Ee)
{
  double tau;
  double tol(1.0e-4);

  setXrayTau(ze,z,Ee,c,this,1);

  //integrate up optical depth
  tau=qromb(getXrayTau,ze,z,tol);

  return tau;
}

//calculate dtau/dzp(zp) for an x-ray observed at z of energy E
double setXrayTau(double z1, double zp, double E1, Cosmology *c1, Astrophysics *a1,int iflag)
{
  static double z,E;
  static Cosmology *c;
  static Astrophysics *a;
  double Ep,nup;
  double dtau;
  double drdz;
  double xe;

  if(iflag==1){
    a=a1;
    c=c1;
    z=z1;
    E=E1;
    return 0.0;
  }

  //NEED A FAST WAY OF EVALUATING THE IONIZATION FRACTION AT ANY Z
  xe=0.0;
  xe=a->getXE(zp);

  Ep=E*(1.0+zp)/(1.0+z);
  nup=Ep/RYDBERG*NULYMANLIMIT;

  drdz=SPEEDOFLIGHT_CGS;
  drdz/=-c->hubbleZ(zp)*UNH*c->getH()*(1.0+zp);

  dtau=a->xsectionPhotoIoniseGen(nup,1,1);       //HI
  dtau+=FHE*a->xsectionPhotoIoniseGen(nup,2,1);  //HeI
  dtau*=drdz*c->nh(zp)*(1.0-xe);

  return dtau;
} 

double getXrayTau(double zp)
{
  return setXrayTau(0.0,zp,0.0,NULL,NULL,0);
}



////////////////////////////////////////////////////////////////////////
//  Calculate Xray heating at a given redshift
///////////////////////////////////////////////////////////////////////

//calculate mean xray heating at redshift z
//Units: xrayflux erg cm^-3 s^-1
double Astrophysics::xrayHeating(double z, double xe)
{
  double xrayheat(0.0);
  double tol(1.0e-4);
  double zstart;
  double zmax(ZSFR);
  //double *result;
  //double tk,xi,xe;
  double z1,z2,dz;
  int i,imax(10);
  
  zmax=min(z+8.0,ZSFR);

  //sum over different levels
  setDJxrayDzDE(20.0,z,0.0,c,this,1);
  getDLambdaXrayDzDE(xe,z,this,1);
  zstart=z;

  //break integral into pieces separated logarithmically
  dz=exp(log(zmax/zstart)/(double)(imax-1));
  xrayheat=0.0;
 
  //  cout<<"beginning loop"<<endl;
  z1=zstart;
  z2=zstart*dz;
  for(i=1;i<=imax;i++){
    xrayheat+=qromb(getDLambdaXrayDz,z1,z2,tol);  
    //    cout<<zstart<<"\t"<<zmax<<"\t"<<z1<<"\t"<<z2<<"\t"<<xrayheat<<endl;
    //cout<<i<<endl;
    z1=z2;
    z2*=dz;
  }

  // cout<<"#"<<zstart<<"\t"<<zmax<<"\t"<<z1<<"\t"<<z2<<"\t"<<xrayheat<<endl;
  
  return xrayheat;
}


//return dLambdadz having integrated over energies
double getDLambdaXrayDz(double zp)
{
  double Emin(RYDBERG);
  double Emax(XrayEmax);
  double dJdz;
  double lEmin,lEmax;
  double tol(1.0e-4);

  lEmin=log(Emin);
  lEmax=log(Emax);

  zpsave=zp;
  dJdz=qromb(getDLambdaXrayDlE,lEmin,lEmax,tol);

  return dJdz;
}

double getDLambdaXrayDlE(double lE)
{
  double E,zp;
  E=exp(lE);
  zp=zpsave;

  return getDLambdaXrayDzDE(E,zp,NULL,0)*E*RYDBERG_CGS/RYDBERG;
}

//returns differential xray heating rate per atom as function of energy E and
//emission redshift zp
//call with zp =z the desired redshift at which the heating occurs
//     and   E = xe the ionization fraction at that z
double getDLambdaXrayDzDE(double E, double zp,Astrophysics *a1, int iflag)
{
  static Astrophysics *a;
  double dJ,dlambda;
  static double fheat(0.0);
  double dlambdaHI,dlambdaHeI;
  double EthHI(RYDBERG);
  double EthHeI(24.59);
  double nu;

  if(iflag==1){
    a=a1;
    fheat=a->fracPrimaryElectron(E,1);
    return 0.0;
  }

  nu=E/RYDBERG*NULYMANLIMIT;
  dJ=getDJxrayDzDE(E,zp);

  //Assumes that ionization fraction of HI and HeI are equal
  //HI contribution
  dlambdaHI=4.0*PI*a->xsectionPhotoIoniseGen(nu,1,1);
  dlambdaHI*=(E-EthHI)*RYDBERG_CGS/RYDBERG;
  //HeI contribution
  dlambdaHeI=FHE*4.0*PI*a->xsectionPhotoIoniseGen(nu,2,1);
  dlambdaHeI*=(E-EthHeI)*RYDBERG_CGS/RYDBERG;
  //HeII
  //assume HeII fraction is negligible

  dlambda=dlambdaHI+dlambdaHeI;
  dlambda*=dJ;
  dlambda*=fheat;

  return dlambda;
}

////////////////////////////////////////////////////////////////////////
//  Calculate \lya flux from xray excitation of HI
///////////////////////////////////////////////////////////////////////

//calculate mean lya flux from xrays at redshift z
//Units: lyafluxXray erg cm^-3 s^-1
double Astrophysics::lyafluxXray(double z, double xe)
{
  double lyaflux(0.0);
  double tol(1.0e-4);
  double zstart;
  double zmax(ZSFR);
  double z1,z2,dz;
  int i,imax(10);
  
  zmax=min(z+8.0,ZSFR);

  //sum over different levels
  setDJxrayDzDE(20.0,z,0.0,c,this,1);
  getDJalphaXrayDzDE(xe,z,this,c,1);
  zstart=z;

  //break integral into pieces separated logarithmically
  dz=exp(log(zmax/zstart)/(double)(imax-1));
  lyaflux=0.0;

  //  cout<<"beginning loop"<<endl;
  z1=zstart;
  z2=zstart*dz;
  for(i=1;i<=imax;i++){
    lyaflux+=qromb(getDJalphaXrayDz,z1,z2,tol);  
    z1=z2;
    z2*=dz;
  }
  
  return lyaflux;
}


//return dJalphaXray having integrated over energies
double getDJalphaXrayDz(double zp)
{
  double Emin(RYDBERG);
  double Emax(XrayEmax);
  double dJdz;
  double lEmin,lEmax;
  double tol(1.0e-4);

  lEmin=log(Emin);
  lEmax=log(Emax);

  zpsave=zp;
  dJdz=qromb(getDJalphaXrayDlE,lEmin,lEmax,tol);

  return dJdz;
}

double getDJalphaXrayDlE(double lE)
{
  double E,zp;
  E=exp(lE);
  zp=zpsave;

  return getDJalphaXrayDzDE(E,zp,NULL,NULL,0)*E*RYDBERG_CGS/RYDBERG;
}

//returns contribution to lya flux from xrays of energy E emitted at redshift 
//zp
//call with zp =z the desired redshift at which the heating occurs
//     and   E = xe the ionization fraction at that z
double getDJalphaXrayDzDE(double E, double zp,Astrophysics *a1, Cosmology *c1,int iflag)
{
  static Astrophysics *a;
  static Cosmology *c;
  static double z;
  double dJ,djalpha;
  static double falpha(0.0);
  double djalphaHI,djalphaHeI(0.0);
  double EthHI(RYDBERG);
  // double EthHeI(24.59);
  double nu;

  if(iflag==1){
    a=a1;
    c=c1;
    z=zp;
    falpha=a->fracPrimaryElectron(E,4);
    return 0.0;
  }

  nu=E/RYDBERG*NULYMANLIMIT;
  dJ=getDJxrayDzDE(E,zp);

  //Assumes that ionization fraction of HI and HeI are equal
  //HI contribution
  djalphaHI=4.0*PI*a->xsectionPhotoIoniseGen(nu,1,1);
  djalphaHI*=(E-EthHI)*RYDBERG_CGS/RYDBERG;
  djalphaHI/=a->hydrogenEnergy(1.0,2.0);
  djalphaHI/=a->hydrogenFreq(1.0,2.0)*c->hubbleZ(z)*UNH*c->getH();

  //HeI contribution
  //  djalphaHeI=FHE*4.0*PI*a->xsectionPhotoIoniseGen(nu,2,1);
  //djalphaHeI*=(E-EthHeI)*RYDBERG_CGS/RYDBERG;
  djalphaHeI=0.0;

  //HeII
  //assume HeII fraction is negligible

  djalpha=djalphaHI+djalphaHeI;
  djalpha*=dJ;
  djalpha*=falpha;
  djalpha*=SPEEDOFLIGHT_CGS/4.0/PI;

  return djalpha;
}





////////////////////////////////////////////////////////////////////////
//  Primary photoelectron fractions
////////////////////////////////////////////////////////////////////////

//Calculate the fraction of a primary electron deposited in
// 1) heating
// 2) ionisation of HI
// 3) ionisation of HeI
// 4) excitation of HI
// 5) excitation of HeI
//Based on Shull and van Steenberg (1985)
//Valid for primary electron energies E0>>100 eV
double Astrophysics::fracPrimaryElectron(double x, int iflag)
{
  double C,a,b,f;

  if(x<0.0) cout<<"Must have positive x in fracPrimaryElectron"<<endl;

  if(iflag==1){
    C=0.9971;
    if(fabs(x-1.0)<1.0e-4) return C;
    a=0.2663;
    b=1.3163;
    f=C*(1.0-pow(1.0-pow(x,a),b));
    return f;
  }

  if(iflag==2){
    C=0.3908;
    a=0.4092;
    b=1.7592;
  }else if(iflag==3){
    C=0.0554;
    a=0.4614;
    b=1.6660;
  }else if(iflag==4){
    C=0.4766;
    a=0.2735;
    b=1.5221;
  }else if(iflag==5){
    C=0.0246;
    a=0.4049;
    b=1.6594;
  }else{
    cout<<"illegal flag in fracPrimaryElectron"<<endl;
    return 0.0;
  }
  if(fabs(x-1.0)<1.0e-4) return 0.0;
  f=C*pow(1.0-pow(x,a),b);
  return f;
}

////////////////////////////////////////////////////////////////
// Star formation rate
///////////////////////////////////////////////////////////////

double Astrophysics::globalSFR(double z)
{
  double sfr;
  int sfr_flag(1);

  if(sfr_flag==1){
    sfr=globalSFRcol(z);
  }else{
    sfr=globalSFRsim(z);
  }    

  return sfr;
}

//Calculate the physical global star formation rate
// Units:  sfr  M_sol Mpc^-3 yr^-1
//
// which mass function I use makes little difference so I use PS
double Astrophysics::globalSFRcol(double z)
{
  double sfr;
  double h,z1;
  h=c->getH();
  z1=z+1.0;

  if(z>ZSFR) return 0.0;

  sfr=dfcdzFast(z);
  sfr*=z1*c->hubbleZ(z)*h*UNH*YEAR;   //dzdt
  sfr*=c->getOmegab()*CRITDENMSOLMPC*h*h*z1*z1*z1;  //rho_b
  sfr*=FSTAR;                             //star formation efficiency

  return sfr;
}

//obtain physical global star formation rate from file
double Astrophysics::globalSFRsim(double z)
{
  static int iflag;
  static Spline sfrSP;
  ifstream fin;
  int i,ndata(91);
  char *file;
  char tag[40];
  double *zz,*sfrIn;
  double datum;

  double sfr;
  double h,z1;
  h=c->getH();
  z1=z+1.0;

  if(iflag==0){
    cout<<"gettting SFR spline data"<<endl;
    iflag++;
    zz=dvector(1,ndata);
    sfrIn=dvector(1,ndata);
    file="simSFR.dat";
    fin.open(file);
    fin.getline(tag,40);
    for(i=1;i<=ndata;i++){
      fin>>zz[i];
      fin>>datum;
      sfrIn[i]=datum;
      fin>>datum;
      sfrIn[i]+=datum;
      cout<<zz[i]<<"\t"<<sfrIn[i]<<endl;
    }
    fin.close();
   
    cout<<"setting SFR spline"<<endl;
    sfrSP.setSplineSP(ndata,zz,sfrIn);

    cout<<"SFR spline set"<<endl;
    cout<<sfrSP.returnValue(25.0)<<endl;

    free_dvector(zz,1,ndata);
    free_dvector(sfrIn,1,ndata);
  }

  if(z>ZSFR) return 0.0;

  sfr=sfrSP.returnValue(z);
  return sfr;
}

////////////////////////////////////////////////////////////////////
//  Compton Heating Rate
////////////////////////////////////////////////////////////////////

// Compton heating for the global temperature evolution
// Heating due to scattering of CMB photons off free electrons
//
// Units: Erg Second^-1 Centimeter^-3
double Astrophysics::heatingCompton(double z, double Tgas, double xfree, Cosmology *c)
{
  double compton;
  double ucmb;
  double ne;
  double Tcmb;

  Tcmb=TCMB*(1.0+z);

  ne=c->nh(z)*xfree;   
  ucmb=RADIATIONA*Tcmb*Tcmb*Tcmb*Tcmb;

  compton=4.0*SIGMATHOMSON/ELECTRONMASS/SPEEDOFLIGHT_CGS;
  compton*=ne*ucmb*BOLTZK*(Tcmb-Tgas);
  return compton;
}

////////////////////////////////////////////////////////////////////
//  Hydrogen recombination rate in early Universe
////////////////////////////////////////////////////////////////////
//From Seager, Sasselov, Scott 1999
//case B recombination rate originally from Hummer (1994)
//Units:  Tgas  K
//        recomb  cm^3 s^-1
double Astrophysics::recombH(double Tgas)
{
   double F(1.14);
   double a(4.309);
   double b(-0.6166);
   double c(0.6703);
   double d(0.5300);
   double t, recomb;

   t=Tgas/1.0e4;

   recomb=F*1e-19;
   recomb*=a*pow(t,b)/(1.0+c*pow(t,d));
   recomb*=1.0e6;  //convert from m^3 to cm^3

   return recomb;
}

////////////////////////////////////////////////////////////////////////
// Spline for the collapse fraction
///////////////////////////////////////////////////////////////////////
// Create a spline of the collapse function for an arbitrary mass function
// initialises spline on first run and runs quickly from then on
// xuse = redshift to evaluate fColl at.
//uses files to save overhead
double Astrophysics::splineFColl(double zuse,Cosmology *c)
{
  static double *ly2;
  static double *lxs, *lys;
  static int nmaxs;
  static int iflag,icount;
  double yp1,ypn;
  double lx0,lxstep,lxmax;
  int n;
  double yout,lxuse;
  int massfcn(MASSFCN);  //0: PS  1: ST  2: Jenkins
  char *file="./fcoll_spline.dat";
  ifstream fin(file);
  ofstream fout;

  lxuse=log(zuse);

  //initialise spline on first pass
  if(iflag==0){
    if(icount>0){
      free_dvector(lxs,1,nmaxs);
      free_dvector(lys,1,nmaxs);
      free_dvector(ly2,1,nmaxs);
    }
    // store information in static variables
    nmaxs=50;
    ly2=dvector(1,nmaxs);
    lxs=dvector(1,nmaxs);
    lys=dvector(1,nmaxs);
    cout <<"initialising fcoll spline"<<endl;

    //if file doesn't exist then initialise by hand
    if(!fin){
      cout<<"calculating fcoll spline"<<endl;
      // range from z=3 to 47.9
      // lx0=1.0;
      lx0=log(ZLOW-0.1);
      lxmax=log(50.0);
      //      lxstep=0.07;
      lxstep=(lxmax-lx0)/(double)(nmaxs-1);
      for(n=1;n<=nmaxs;n++){
	lxs[n]=lx0+lxstep*(double)(n-1);
	lys[n]=log(c->fColl(exp(lxs[n]),coolMass(c,exp(lxs[n])),massfcn));
	//	cout<<n<<"\t"<<exp(lxs[n])<<"\t"<<lys[n]<<endl;
      }
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
    cout<<"fcoll spline initialised"<<endl;
    return exp(yout);
  }

  // apply cubic spline
  if(iflag==1){
    if((lxuse<lxs[1])||(lxuse>lxs[nmaxs])) {
     cout <<"lx exceeds limits in splineFColl: "<<exp(lxuse)<<endl;
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

  cout << "Error in splineFColl"<<endl;
  return 0.0;
}

//Calculate dfColl/dz using spline of fColl
double Astrophysics::dfcdzFast(double z)
{
  double dfcdz;

  dfcdz=fabs((splineFColl(z+0.001,c)-splineFColl(z-0.001,c))/0.002);

  return dfcdz;
}

///////////////////////////////////////////////////////////////////
// Spline the free electron fraction history
///////////////////////////////////////////////////////////////////

//create and use a spline of the Window function
double Astrophysics::splineXFree(double zuse, double z[], double xe[], int nmax, int iflag)
{
  static double *y2;
  static double *xs, *ys;
  static int nmaxs;
  int n;
  double yout;
  double yp1(1.0e30);
  double ypn(1.0e30);
  double xuse;
  
  xuse=zuse;

  //initialise spline on first pass
  if(iflag==0){
    // store information in static variables
    nmaxs=nmax;
    y2=dvector(1,nmax);
    xs=dvector(1,nmax);
    ys=dvector(1,nmax);
    //  for(n=1;n<=nmax;n++){
    for(n=1;n<=nmax;n++){
      xs[n]=z[nmax-n+1];
      ys[n]=log(xe[nmax-n+1]);
      //    cout<<xs[n]<<"\t"<<ys[n]<<endl;
    }
    //initialise spline
    spline(xs,ys,nmaxs,yp1,ypn,y2);
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    return exp(yout);
  }
  // apply cubic spline
  if(iflag==1){
    if((xuse<min(xs[1],xs[nmaxs]))||(xuse>max(xs[1],xs[nmaxs]))) {
      cout <<"x exceeds limits in splineXfree"<<endl;
      cout <<xuse<<"\t"<<xs[1]<<"\t"<<xs[nmaxs]<<endl;
      return 0.0;
    }
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    return exp(yout);
  }

  if(iflag==2){
    free_dvector(xs,1,nmaxs);
    free_dvector(ys,1,nmaxs);
    free_dvector(y2,1,nmaxs);
    return 0.0;
  }

  cout << "Error in splineXfree"<<endl;
  return 0.0;
}

// easy to use function based on the proceeding spline
double Astrophysics::getXFree(double zuse)
{
  return splineXFree(zuse,NULL,NULL,0,1);
}

////////////////////////////////////////////////////////////////
// Integrate Thermal and Ionization history
////////////////////////////////////////////////////////////////

// Integrate up the history of T_IGM and x_e
// Output full history in "Thistory_save.dat"
// Return values of T and x_e at zin
//
//Units: T  K
void Astrophysics::globalHistory(void)
{
  double z,zstart(ZSFR);
  double zstep(0.25), zmin(ZLOW);
  double *y;
  double *dy;
  double globalxe,globalT,globalxi;
  int nvar(3);
  int nok,nbad;
  double hstep(1.0);
  double hmin(0.001);
  double eps(1e-4);
  double z0,tk0,xi0;
  ofstream fout;
  char *file;

  double *zs, *tks, *xes, *xis;
  double datastore[500][4];
  int i,ndata(1);

  y=dvector(1,nvar);
  dy=dvector(1,nvar);
   
  setDerivsHistory(zmin,y,dy,c,this,1);

  cout<<"calculating global history"<<endl;

  //call RECFAST and extract T_IGM and x_e history before reionization

  //use output of recfast to set initial conditions
  useRECFAST(zstart,y);
  z0=y[1];
  tk0=y[2];
  xi0=y[3];

  zstart=z0;
  y[1]=tk0;
  //y[2]=xi0;
  y[2]=0.0;  //initially there are no bubbles
  y[3]=xi0;

  //store initial conditions then loop over decreasing redshift
  datastore[ndata][0]=zstart;
  datastore[ndata][1]=y[1];
  datastore[ndata][2]=y[2];
  datastore[ndata][3]=y[3];
  while(zstart>zmin && zstart>0.0){
    ndata++;
    z=zstart-zstep;

    odeint(y,nvar,zstart,z,eps,hstep,hmin,&nok,&nbad,derivsHistory,rkqs);

    datastore[ndata][0]=z;
    datastore[ndata][1]=y[1];
    datastore[ndata][2]=y[2];
    datastore[ndata][3]=y[3];
    //    cout<<z<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<y[3]<<endl;
    zstart=z;
  }

  zs=dvector(1,ndata);
  tks=dvector(1,ndata);
  xis=dvector(1,ndata);
  xes=dvector(1,ndata);

  // output calculated information in ascending redshift order
  file="Thistory_save.dat";
  fout.open(file);
  for(i=ndata;i>0;i--){
    fout<<datastore[i][0]<<"\t"<<datastore[i][1]<<"\t"<<datastore[i][2]<<"\t"<<datastore[i][3]<<endl;

    zs[ndata-i+1]=datastore[i][0];
    tks[ndata-i+1]=datastore[i][1];
    xis[ndata-i+1]=datastore[i][2];
    xes[ndata-i+1]=datastore[i][3];
  }
  fout.close();

  //set splines for future use
  //  cout<<ndata<<endl;
  globalTK.setSplineSP(ndata,zs,tks);
  globalXI.setSplineSP(ndata,zs,xis);
  globalXE.setSplineSP(ndata,zs,xes);

  //tidy memory
  free_dvector(y,1,nvar);  
  free_dvector(dy,1,nvar);
  free_dvector(zs,1,ndata);
  free_dvector(tks,1,ndata);
  free_dvector(xis,1,ndata);
  free_dvector(xes,1,ndata);

  cout<<"global history calculated successfully"<<endl;
  globalHistoryFlag=1;
  ZREION=findZReion();
  cout<<"reionization at z="<<ZREION<<endl;
}

// Return values of T and x_e at zin
// uses full integration, but doesn't work out full history
//
//Units: T  K
void Astrophysics::getTIGMnum(double zin, double result[])
{
  double zstart(ZSFR);

  double *y;
  double *dy;
  double globalxe,globalT,globalxi;
  int nvar(3);
  int nok,nbad;
  double hstep(1.0);
  double hmin(0.001);
  double eps(1e-4);
  double z0,tk0,xi0;
  
  y=dvector(1,nvar);
  dy=dvector(1,nvar);
   
  setDerivsHistory(zin,y,dy,c,this,1);

  //If zin >SFR extract tk,xi,xe from RECFAST history
  if(zin>zstart){
    getRECFAST(zin,y);
    result[1]=y[2];
    //    result[2]=y[3];
    result[2]=0.0;  //no HII regions before star formation begins
    result[3]=y[3];
    return;
  }

  //use output of recfast to set initial conditions
  useRECFAST(zstart,y);
  z0=y[1];
  tk0=y[2];
  xi0=y[3];

  zstart=z0;
  y[1]=tk0;
  //y[2]=xi0;
  y[2]=0.0;  //initially there are no bubbles
  y[3]=xi0;
  odeint(y,nvar,zstart,zin,eps,hstep,hmin,&nok,&nbad,derivsHistory,rkqs);

  globalT=y[1];
  globalxi=y[2];
  globalxe=y[3];
  
  free_dvector(y,1,nvar);  
  free_dvector(dy,1,nvar);

  result[1]=globalT;
  result[2]=globalxi;
  result[3]=globalxe;
}

// From global history find out when the Universe became ionized
//
double Astrophysics::findZReion()
{
  ifstream fin;
  char *file;
  double z, tk, xi, xe;
  double zmark;

  file="Thistory_save.dat";
  fin.open(file);
  while(!fin.eof()){
    fin>>z>>tk>>xi>>xe;
    //    cout<<z<<"\t"<<tk<<"\t"<<xi<<"\t"<<xe<<endl;
    if(fabs(xi-1.0)<1.0e-4) zmark=z;
  }

  return zmark;
}

//More detailed finding of reionization redshift
double Astrophysics::findZReionDetail()
{
	double zreion,z;
	double global_step(0.5);  //step size in globalHistory()
	double step(0.01);
	double xi;
	double *result;

	result=dvector(1,3);

	zreion=findZReion();
	cout<<"zreion="<<zreion<<endl;
	z=zreion+global_step+step;
	while(z>=zreion){
		z-=step;
		getTIGMnum(z,result);
		xi=result[2];
		cout<<z<<"\t"<<xi<<endl;
		if(fabs(xi-1.0)<1.0e-4) break;
	}

	free_dvector(result,1,3);
	return z;
}


// Get IGM evolution using splines of global history
// with simple handling of the reionized period
//
void Astrophysics::getTIGM(double zin, double result[])
{
  if(globalHistoryFlag==0) globalHistory();

  if(zin>ZSFR){  //before star formation begins
    getRECFAST(zin,result);
    result[1]=result[2];
    result[2]=0.0;  //no HII regions before star formation begins
  }else if(zin<ZREION){  //after reionization completes
    result[1]=3.0e4;  //photo-heating dominates
    result[2]=1.0;    //fully ionized
    result[3]=1.0;
  }else{            //before reionization but after star formation begins
    result[1]=globalTK.returnValue(zin);
    result[2]=globalXI.returnValue(zin);
    result[3]=globalXE.returnValue(zin);
    if(result[2]>1.0) result[2]=1.0;
  }
}

// Get IGM evolution using splines of global history
//
//
void Astrophysics::getTIGMcont(double zin, double result[])
{
  if(zin>ZSFR){  //before star formation begins
    getRECFAST(zin,result);
    result[1]=result[2];
    result[2]=0.0;  //no HII regions before star formation begins
  }else{            //before reionization but after star formation begins
    result[1]=globalTK.returnValue(zin);
    result[2]=globalXI.returnValue(zin);
    result[3]=globalXE.returnValue(zin);
    if(result[2]>1.0) result[2]=1.0;
  }
}

// function to calculate derivatives for globalHistory
void setDerivsHistory(double z, double y[],double dy[],Cosmology *c1,Astrophysics *a1,int iflag)
{
  static Cosmology *c;
  static Astrophysics *a;
  double dtg,zz,tg,dx;
  double heat;
  double CT,dzdt,ntot;
  double xi,xfree;
  double zeta;  // ionisation efficiency parameter
  double ion,recomb;
  double alphaA(4.2e-13);  //case-A recombination: cm^3 s^-1 (T=10^4K)
  double C(0.0);  //Clumping factor - NEED MODEL FOR THIS!!!!!!
                  //constant C fails badly at high x_e!!


  if(iflag==1){
    c=c1;
    a=a1;
  }else{

    zeta=a->getZeta();
    C=CLUMPING;

    tg=y[1];
    if(y[2]>1.0) y[2]=1.0;
    xi=y[2];
    if(y[3]>1.0) y[3]=1.0;
    xfree=y[3];

    if(xi>1.0e-4){
      a->initPVParamMEHR(z);
      C=a->clumpingMEHR(z,xi);
      //      cout<<z<<"\t"<<xi<<"\t"<<C<<endl;
      if(z<4.0) C=20.0;  //hack to avoid C problems at low z
    }
    
  ///////////////////////////
  //first handle the temperature evolution

    zz = z+1.0;
    dtg = 2.0*tg/zz;  //adiabatic expansion
    
    //Add in heating terms
    ntot=c->nh(z)*(1.0+xfree+FHE);  //n=n_H+n_e+n_He
    CT=1.5*BOLTZK*ntot;  //specific heat capacity in ergs K^{-1} cm^{-3}
    dzdt=zz*c->hubbleZ(z)*UNH*c->getH();
    
    heat=0.0;
    heat=a->heatingCompton(z,tg,xfree,c);  //compton 
    heat+=a->lumTotStarBurst(a->globalSFR(z))/MPC/MPC/MPC*a->fracPrimaryElectron(xfree,1); //xrays
    //   heat+=a->xrayHeating(z,xfree); //xrays
    
    heat/=dzdt*CT; 
    dtg-= heat;
    dy[1]= dtg;
    if(tg>1.0e5) dy[1]=0.0;  //clamp max temp to 100,000K
  
    //////////////////////////////
    //Second handle ionization evolution, x_i
    ion=zeta*a->dfcdzFast(z); 
    ion*=(1.0-xfree);       //Account for free electron fraction
    recomb=alphaA*C*xi*xi*c->nh(z);  
    recomb/=dzdt;
    
    dx=ion-recomb;
    
    dy[2]= -dx;
    if(fabs(xi-1.0)<1.0e-4) dy[2]=0.0; //cheat in case that xfree>1 
    

    //////////////////////////////
    //Third handle ionization evolution, x_e in IGM
    //ionizing radiation is all from x-rays
    ion=a->lumTotStarBurst(a->globalSFR(z))/MPC/MPC/MPC*a->fracPrimaryElectron(xfree,2)/RYDBERG_CGS; //xrays
    //     ion=a->xrayHeating(z,xfree)*a->fracPrimaryElectron(xfree,2)/a->fracPrimaryElectron(xfree,1)/RYDBERG_CGS; //xrays

    ion/=c->nh(z);   //ionizing photons per hydrogen atom
    ion*=(1.0-xfree); //Account for free electron fraction
    ion/=dzdt;
    
    //recombination in the bulk is temperature dependent and case B
    recomb=a->recombH(tg)*xfree*xfree*c->nh(z);
    recomb/=dzdt;
    
    dx=ion-recomb;
    
    dy[3]= -dx;
    if(fabs(xi-1.0)<1.0e-4) dy[3]=0.0; //cheat in case that xi>1, so after reionization
    
  }
}

//dummy function to get derivative for odeint for globalHistory
void derivsHistory(double z, double y[], double dy[])
{
  setDerivsHistory(z,y,dy,NULL,NULL,0);
}
////////////////////////////////////////////////////////////////
// Get Tk, xi and xe quickly
////////////////////////////////////////////////////////////////
//should change these to calculate history once and then get
//values from a spline
double Astrophysics::getTK(double z)
{
  double tk;
  double *result;
  result=dvector(1,3);
  getTIGM(z,result);
  tk=result[1];
  free_dvector(result,1,3);
  return tk;
}

double Astrophysics::getXI(double z)
{
  double xi;
  double *result;
  result=dvector(1,3);
  getTIGM(z,result);
  xi=result[2];
  free_dvector(result,1,3);
  return xi;
}

double Astrophysics::getXE(double z)
{
  double xe;
  double *result;
  result=dvector(1,3);
  getTIGM(z,result);
  xe=result[3];
  free_dvector(result,1,3);
  return xe;
}


////////////////////////////////////////////////////////////////
// Calculate tau for a specified ionization history
////////////////////////////////////////////////////////////////


double Astrophysics::getTauCMB()
{
  double tau;
  double zmin(0.01),zmax(40);
  double tol(1.0e-4);

  setTauInt(zmin,ZREION,c,this,0);

  tau=qromb(getTauInt,zmin,zmax,tol);

  return tau;
}

//integrand for tau integral
double setTauInt(double z, double zri1,Cosmology *c1, Astrophysics *a1,int iflag)
{
  double dtaudz;
  static Cosmology *c;
  static Astrophysics *a;
  static double zri;
  double zHe(3.0);
  double xi, xe, xmean;

  if(iflag==0){
    zri=zri1;
    c=c1;
    a=a1;
    return 0.0;
  }

  dtaudz=c->nh(z)*SPEEDOFLIGHT_CGS*SIGMATHOMSON;
  dtaudz/=(1.0+z)*c->hubbleZ(z)*UNH*c->getH();
  if(z>zri){
    xi=a->getXI(z);
    xe=a->getXE(z);
    xmean=xi+(1.0-xi)*xe;  //volume weighted ionized fraction
    dtaudz*=xmean; 
  }
  if(z<zHe) dtaudz*=1.0+2.0*FHE;  //account for Helium reionization

  return dtaudz;
} 

//dummy function for tau integral
double getTauInt(double z)
{
  return setTauInt(z,0.0,NULL,NULL,1);
}

////////////////////////////////////////////////////////////////////////
//  Calculate the clumping factor C
////////////////////////////////////////////////////////////////////////

// Model taken from Miralde-Escude, Haehnelt, Rees (2000)

//  calcualte clumping factor at redshift z given bubble filling fraction xi
//
double Astrophysics::clumpingMEHR(double z, double xi)
{
  double clump;
  double ldeltai;
  double ldeltaMin(-14.0);
  double tol(1.0e-4);

  if(xi<1.0e-4) return 0.0;
  if(fabs(xi-1.0)<1.0e-8) return 1.0e1;  
  // need to handle things carefully

  ldeltai=log(deltaIMEHR(z,xi));
  //  clump=qromb(getlDelta3PVMEHR,ldeltaMin,ldeltai,tol);
  clump=clumpingIntMEHR(ldeltai);

  return clump;
}

// calculate integral of P_V \delta^2
//
double Astrophysics::getClumpingMEHR(double Deltai, double z)
{
  double clump;
  double lDeltai;
  double lDeltaMin(-14.0);
  double tol(1.0e-4);

  lDeltai=log(Deltai);
  //Set up parameters for P_V distribution
  //  initPVParamMEHR(z);
  // Calculate the ionization threshold Delta_i
  //  Deltai=2.0;
  //calculate the clumping factor
  clump=qromb(getlDelta3PVMEHR,lDeltaMin,lDeltai,tol);
  //cout<<"clump:"<<clump<<"\t"<<clumpingIntMEHR(lDeltai)<<endl;
  clump=clumpingIntMEHR(lDeltai);
  return clump;
}

// calculate integral of P_V \delta^2
//
double Astrophysics::getVolumeMEHR(double Deltai, double z)
{
  double fv;
  double lDeltai;
  double lDeltaMin(-14.0);
  double tol(1.0e-4);

  lDeltai=log(Deltai);
  //Set up parameters for P_V distribution
  //  initPVParamMEHR(z);
  // fv=qromb(getlDeltaPVMEHR,lDeltaMin,lDeltai,tol);
  fv=volumeIntMEHR(lDeltai);

  return fv;
}

// Calculate critical overdensity corresponding to F_V(\Delta_i)=x_i
// This can then be used to get clumping factor
//
double Astrophysics::deltaIMEHR(double z, double xi)
{
  double tol(1.0e-4);
  double ldelta1(-14.0), ldelta2(14-1.0e-4);
  double ldeltai;
  
  //set up parameters for P_V distribution
  //initPVParamMEHR(z);
  //root find to get the value of delta_i
  setDeltaIMEHR(z,xi,this,1);
  ldeltai=zriddrSimp(dummyDeltaIMEHR,ldelta1,ldelta2,tol);

  return exp(ldeltai);
}

double setDeltaIMEHR(double ldeltai, double xi1, Astrophysics *a1,int iflag)
{
  static double xi;
  static Astrophysics *a;
  double fv;
  double tol(1.0e-4);
  double ldeltaMin(-14.0);

  if(iflag==1){
    xi=xi1;
    a=a1;
    return 0.0;
  }
  
  //  cout<<ldeltai<<endl;
  if(ldeltai<ldeltaMin) return -xi;
  //  fv=qromb(getlDeltaPVMEHR,ldeltaMin,ldeltai,tol);
  fv=a->volumeIntMEHR(ldeltai);
  // cout<<"vol:"<<fv<<"\t"<<a->volumeIntMEHR(ldeltai)<<endl;
  //cout<<fv<<"\t"<<xi<<endl;
  return fv-xi;
}

double dummyDeltaIMEHR(double ldeltai)
{
  return setDeltaIMEHR(ldeltai,0.0,NULL,0);
}

void Astrophysics::initPVParamMEHR(double z)
{
  static double zcurrent;
  double x1(0.1), x2(2.0);
  double tol(1.0e-4);
  double zz;
  double A(1.0),C0(0.1),beta,delta0;
  double lDeltaMax(14.0);
  double lDeltaMin(-14.0);

  int i,nspline(80);
	nspline=160;
  double ldstep;
  double *ldvec, *fvvec, *clumpvec;

  zz=z+1.0;

  if(fabs(z-zcurrent)<1.0e-4) return;  //already initiated at correct z
  zcurrent=z;

  delta0=7.61/zz;

  beta=2.23;
  if(z>2.99) beta=2.35;
  if(z>3.99) beta=2.48;
  if(z>5.99) beta=2.50;

  setPVNormMEHR(beta,C0,delta0,1);

  C0=zriddrSimp(getPVNormMEHR,x1,x2,tol);
 
  setPVMEHR(1.0,1.0,beta,C0,delta0,1);
  //  cout<<C0<<"\t"<<qromb(getlDeltaPVMEHR,lDeltaMin,lDeltaMax,tol)<<"\t"<<qromb(getlDelta2PVMEHR,lDeltaMin,lDeltaMax,tol)<<endl;
  A=1.0/qromb(getlDeltaPVMEHR,lDeltaMin,lDeltaMax,tol);

  setPVMEHR(1.0,A,beta,C0,delta0,1);

  //  cout<<z<<"\t"<<A<<"\t"<<delta0<<"\t"<<beta<<"\t"<<C0<<endl;
  //now initialise splines for volume and clumping
  ldvec=dvector(1,nspline);
  fvvec=dvector(1,nspline);
  clumpvec=dvector(1,nspline);

  // splines are annoyingly bad (1.0e-3 error) at low delta<-0.4
  // acceptable at values larger than that
  //
  ldstep=(lDeltaMax-lDeltaMin)/(double)(nspline-1);
  for(i=1;i<=nspline;i++){
    ldvec[i]=lDeltaMin+(double)(i-1)*ldstep;
    fvvec[i]=log(qromb(getlDeltaPVMEHR,lDeltaMin,ldvec[i],tol)+1.0e-300);
    clumpvec[i]=log(qromb(getlDelta3PVMEHR,lDeltaMin,ldvec[i],tol)+1.0e-300);
    //    cout<<ldvec[i]<<"\t"<<fvvec[i]<<"\t"<<clumpvec[i]<<endl;
  }
  
  MEHRclumping.setSplineSP(nspline,ldvec,clumpvec);
  MEHRvolume.setSplineSP(nspline,ldvec,fvvec);

  free_dvector(ldvec,1,nspline);
  free_dvector(fvvec,1,nspline);
  free_dvector(clumpvec,1,nspline);
}

double Astrophysics::getMEHRindex(double z)
{
  double beta;
  beta=2.23;
  if(z>2.99) beta=2.35;
  if(z>3.99) beta=2.48;
  if(z>5.99) beta=2.50;
  return beta;
}

//Dummy function for root finding in clumping routines
double getPVNormMEHR(double C0)
{
  return setPVNormMEHR(0.0,C0,0.0,0);
}


//Function for root finding algorithm in clumping to get C0
double setPVNormMEHR(double beta1, double C0, double delta01,int iflag)
{
  double tol(1.0e-4);
  double lDeltaMax(-14.0);
  double lDeltaMin(12.0);
  double err;
  static  double beta,delta0;

  if(iflag==1){
    beta=beta1;
    delta0=delta01;
    setPVMEHR(0.0,1.0,beta,C0,delta0,1);
    return 0.0;
  }

  setPVMEHR(0.0,1.0,beta,C0,delta0,1);
  //when difference between \int P_V and \int Delta P_V is zero we've found C0
  err=qromb(getlDeltaPVMEHR,lDeltaMin,lDeltaMax,tol)-qromb(getlDelta2PVMEHR,lDeltaMin,lDeltaMax,tol);
  return err;
}

//probability distribution for \Delta 
double PVMEHR(double Delta, double A, double beta, double C0, double delta0)
{
  double temp;

  temp=pow(Delta,-2.0/3.0)-C0;
  temp*=-temp;
  temp/=2.0*(2.0*delta0/3.0)*(2.0*delta0/3.0);
  temp=exp(temp);
  temp*=pow(Delta,-beta);
  temp*=A;

  return temp;
}

//set dummy function for PV in clumping function
double setPVMEHR(double Delta, double A1, double beta1, double C01, double delta01, int iflag)
{
  double pv;
  static double A,beta,C0,delta0;

  if(iflag==1){
    A=A1;
    beta=beta1;
    C0=C01;
    delta0=delta01;
    return 0.0;
  }

  pv=PVMEHR(Delta,A,beta,C0,delta0);
  
  return pv;
}

//dummy function for PV in clumping function
double getPVMEHR(double Delta)
{
  return setPVMEHR(Delta,0.0,0.0,0.0,0.0,0);
}

//dummy function for Delta*PV in clumping function
double getDeltaPVMEHR(double Delta)
{
  return Delta*setPVMEHR(Delta,0.0,0.0,0.0,0.0,0);
}

//dummy function for Delta*PV in clumping function
double getDelta2PVMEHR(double Delta)
{
  return Delta*Delta*setPVMEHR(Delta,0.0,0.0,0.0,0.0,0);
}

//dummy function for Delta*PV in clumping function
double getlDelta3PVMEHR(double lDelta)
{
  double Delta;
  Delta=exp(lDelta);
  return Delta*Delta*Delta*setPVMEHR(Delta,0.0,0.0,0.0,0.0,0);
}

//dummy function for Delta*PV in clumping function
double getlPVMEHR(double lDelta)
{
  double Delta;
  Delta=exp(lDelta);
  return setPVMEHR(Delta,0.0,0.0,0.0,0.0,0);
}

//dummy function for Delta*PV in clumping function
double getlDeltaPVMEHR(double lDelta)
{
  double Delta;
  Delta=exp(lDelta);
  return Delta*setPVMEHR(Delta,0.0,0.0,0.0,0.0,0);
}

//dummy function for Delta*PV in clumping function
double getlDelta2PVMEHR(double lDelta)
{
  double Delta;
  Delta=exp(lDelta);
  return Delta*Delta*setPVMEHR(Delta,0.0,0.0,0.0,0.0,0);
}


///////////////////////////////////////////////////////////////
//  RECFAST calling program
//////////////////////////////////////////////////////////////
// Following use modified version of RECFAST which outputs
// z  x_e  tmat  trad


void Astrophysics::callRECFAST(void)
{
  char *file="inRECFAST.dat";
  ofstream fout;
  system("rm ./outRECFAST.dat");
  cout<<"calling RECFAST"<<endl;
  fout.open(file);
  fout<<"outRECFAST.dat"<<endl;
  fout<<c->getOmegab()<<endl;
  fout<<c->getOmega0()<<endl;
  fout<<c->getLambda0()<<endl;
  fout<<c->getH()*100.0<<endl;
  fout<<TCMB<<endl;
  fout<<YP<<endl;
  fout.close();

  system("./recfast <inRECFAST.dat");

  cout<<"RECFAST called successfully"<<endl;
}

//Find values in RECFAST output which are closest to and larger than
// the redshift z
void Astrophysics::useRECFAST(double z, double result[])
{
  static int iflag;
  ifstream fin;
  double zcurrent(1.0e5);
  double xe,tmat,trad;
  double zold,xeold,tmatold,tradold;
  char *file="outRECFAST.dat";

  if(iflag==0) callRECFAST();
  iflag++;

  fin.open(file);
  while(zcurrent>z){
    zold=zcurrent;
    xeold=xe;
    tmatold=tmat;
    trad=tradold;
    fin>>zcurrent>>xe>>tmat>>trad;
  }
  fin.close();
  //  cout<<z<<"\t"<<zcurrent<<"\t"<<xe<<"\t"<<tmat<<"\t"<<trad<<endl;

  //  result[1]=zcurrent;
  //result[2]=tmat;
  //result[3]=xe;
  result[1]=zold;
  result[2]=tmatold;
  result[3]=xeold;

}

//Find values for Tk, xi,xe from splined RECFAST data
void Astrophysics::getRECFAST(double z, double result[])
{
  double tk,xi;

  tk=splineRECFASTT(z);
  xi=splineRECFASTX(z);

  result[1]=z;
  result[2]=tk;
  result[3]=xi;

}

//Create and return splined values of T from RECFAST history
double Astrophysics::splineRECFASTT(double zuse)
{
  static double *y2;
  static double *xs, *ys;
  static int nmaxs;
  static int icount;
  static int iflag;
  static int useflag;
  int n,nmax(1000);
  double yout;
  double yp1(1.0e30);
  double ypn(1.0e30);
  double xuse;
  char *file="./splineRECFASTT.dat";
  char *fileREC="./outRECFAST.dat";
  ifstream fin(file);
  ifstream finREC;
  ofstream fout;
  double tk,xi,trad;
  
  if(useflag==0){ 
    callRECFAST();
    useflag++;
  }

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
    cout <<"initialisaing RECFASTT spline"<<endl;
    //if can't initiate from file so get from RECFAST file
    if(!fin){
      //////////
      double z;
      finREC.open(fileREC);
      cout<<"calculating RECFASTT spline"<<endl;
      //    for(n=1;n<=nmaxs;n++){
      for(n=1000;n>=1;n--){
	finREC>>z>>xi>>tk>>trad;
	xs[n]=z;
	ys[n]=tk;
      }
      finREC.close();

      spline(xs,ys,nmaxs,yp1,ypn,y2);
      fout.open(file);
      for(n=1;n<=nmaxs;n++){
	fout << xs[n]<<"\t"<<ys[n]<<"\t"<<y2[n]<<endl;
      }
      fout.close();

      ///////////
    }else{
      for(n=1;n<=nmaxs;n++){
	fin >>xs[n] >>ys[n] >>y2[n];
      }
      fin.close();
    }
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    icount++;
    iflag++;
    cout <<"RECFASTT spline initiated"<<endl;
    return yout;
  }
  // apply cubic spline
  if(iflag==1){
    if((xuse<xs[1])||(xuse>xs[nmaxs])) {
      cout <<"x exceeds limits in RECFASTT spline"<<endl;
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

  cout << "Error in RECFASTT spline"<<endl;
  return 0.0;
}

//Create and return splined values of X_e from RECFAST history
double Astrophysics::splineRECFASTX(double zuse)
{
  static double *y2;
  static double *xs, *ys;
  static int nmaxs;
  static int icount;
  static int iflag;
  static int useflag;
  int n,nmax(1000);
  double yout;
  double yp1(1.0e30);
  double ypn(1.0e30);
  double xuse;
  char *file="./splineRECFASTX.dat";
  char *fileREC="./outRECFAST.dat";
  ifstream fin(file);
  ifstream finREC;
  ofstream fout;
  double tk,xi,trad;
  
  if(useflag==0){ 
    callRECFAST();
    useflag++;
  }

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
    cout <<"initialisaing RECFASTX spline"<<endl;
    //if can't initiate from file so get from RECFAST file
    if(!fin){
      //////////
      double z;
      finREC.open(fileREC);
      cout<<"calculating RECFASTX spline"<<endl;
      for(n=1000;n>=1;n--){
	finREC>>z>>xi>>tk>>trad;
	xs[n]=z;
	ys[n]=xi;
      }
      finREC.close();

      spline(xs,ys,nmaxs,yp1,ypn,y2);
      fout.open(file);
      for(n=1;n<=nmaxs;n++){
	fout << xs[n]<<"\t"<<ys[n]<<"\t"<<y2[n]<<endl;
      }
      fout.close();

      ///////////
    }else{
      for(n=1;n<=nmaxs;n++){
	fin >>xs[n] >>ys[n] >>y2[n];
      }
      fin.close();
    }
    splint(xs,ys,y2,nmaxs,xuse,&yout);
    icount++;
    iflag++;
    cout <<"RECFASTX spline initiated"<<endl;
    return yout;
  }
  // apply cubic spline
  if(iflag==1){
    if((xuse<xs[1])||(xuse>xs[nmaxs])) {
      cout <<"x exceeds limits in RECFASTX spline"<<endl;
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

  cout << "Error in RECFASTX spline"<<endl;
  return 0.0;
}

///////////////////////////////////////////////////////////////////////
// Halo cutoff functions
///////////////////////////////////////////////////////////////////////
/* Jeans mass; BL01 eq. 41 */
double Astrophysics::jeansMassFull(double z, double tk)
{
  double cs,mj,tff,z3;
  //assume neutral gas for molecular weight of gas

  z3=pow(1.0+z,3.0);
  tff=sqrt(4.0*PI*NEWTONCONSTANT*c->getOm0hh()*CRITDENSITY*z3);
  cs=sqrt(5.0/3.0*BOLTZK*tk/(MUB*PROTONMASS));
  mj=pow(PI*cs/tff,3.0);
  mj*=4.0*PI/3.0*CRITDENSITY*c->getOm0hh()*z3/SOLARMASS;

  return mj;
}

// Filtered mass, accounting for full thermal history.  See Gnedin (2000) 
// only valid during matter dominated regime
double Astrophysics::filterMassFull(double z)
{
  double zt;
  double temp,mass;
  double zstart(1000.0);
  double tol(1.0e-4);
  double step(0.05);

  setFilterMassKernel(z,this,1);
  if(z<ZREION && fabs(z-ZREION)>step){
    //tk discontinuous at ZREION, so split integral in two
    temp=qromb(filterMassKernel,ZREION+step,zstart,tol); 
    temp+=step*(filterMassKernel(ZREION+step)+filterMassKernel(ZREION-step))/2.0;
    temp+=qromb(filterMassKernel,z,ZREION-step,tol);
  }else{
    temp=qromb(filterMassKernel,z,zstart,tol);
  }
  mass=pow(temp,1.5);

  return mass;
}

// dummy function for filter mass calculation
double setFilterMassKernel(double z, Astrophysics *a1, int iflag)
{
  static Astrophysics *a;
  static double zuse;
  double mf,mj,tk;
  double zz;

  if(iflag==1){
    zuse=z;
    a=a1;
    return 0.0;
  }

  zz=z+1.0;
  tk=a->getTK(z);
  mj=a->jeansMassFull(z,tk);
  mf=3.0*(zuse+1.0)*pow(mj,2.0/3.0)/zz/zz*(1.0-sqrt((zuse+1.0)/zz));
  return mf;
}

double filterMassKernel(double z)
{
  return setFilterMassKernel(z,NULL,0);
}


////////////////////////////////////////////////////////////////////////
//  Calculate 21 cm Radio flux at a given redshift
///////////////////////////////////////////////////////////////////////
//calculate mean comoving Radio flux at redshift z
// Reference: Barkana and Loeb (2005) detecting
//Units: radioflux cm^2 s^-1 Hz^-1 ster^-1
//       nu  Hz
double Astrophysics::radioFlux(double nu, double z)
{
  double radioflux;
  double tol(1.0e-4);
  double zstart,zmax;
  setDJradioDz(nu,z,z,c,this,1);

  if(z>ZSFR) return 0.0;

  zstart=z;
  zmax=ZSFR;
  radioflux=qromb(getDJradioDz,zstart,zmax,tol);
  
  return radioflux;
}

// Function to calculate the radioEmission, eta
// Units: nu s^-1
//        sE cm^-3 s^-1 Hz^-1
// might be better to make this per unit SFR
double Astrophysics::radioEmission(double nu, double z)
{
  double emission0(2.0e-4);
  double spec_alpha(-0.8);
  double zz;
  double sE(0.0);
  double sfr;
  zz=1.0+z;
  
  // synch emission
  //sE=emission0;
  //sE *= pow(nu/1.0e9,spec_alpha-1.0);

  //Bressan et al updated emissivity
  emission0=5.80e-5;
  sE=0.06/1.44*pow(nu/1.49e9,-1.5);
  sE+=1.38/1.44*pow(nu/1.49e9,spec_alpha-1.0);
  sE*=emission0;

  //free-free from ISM 
  sE+=1.86464e-9*(1-FESC)*NION*pow(nu/1.0e9,-1.0);
  //free-free from IGM - suppressed by FESC/(1-FESC) over ISM contribution
  sE+=1.86464e-9*FESC*NION*pow(nu/1.0e9,-1.0);

  sfr= globalSFR(z)/zz/zz/zz/PROTONMASS;  //comoving SFR rate
  sfr/=1.0e9*KPC*KPC*KPC*YEAR/SOLARMASS;

  sE*=sfr;

  //radio loud AGN
  //  cout<<"sE="<<nu<<"\t"<<sE<<"\t";
  sE+=emissivityRadioAGN(nu,z)/(nu*PLANCKCONSTANT)/MPC/MPC/MPC;
  //cout<<emissivityRadioAGN(nu,z)/(nu*PLANCKCONSTANT)/MPC/MPC/MPC<<endl;

   return sE;
}

// Function to calculate the differential radio flux on a gas element at 
// redshift z, from a source at zp
double setDJradioDz(double nu1, double z1, double zp, Cosmology *c1, Astrophysics *a1,int iflag)
{
  static double nu0;
  static double z;
  static Cosmology *c;
  static Astrophysics *a;
  double dJdz;
  double nup;

  if(iflag==1){
    nu0=nu1;
    z=z1;
    a=a1;
    c=c1;
  }
  if(zp>ZSFR) return 0.0;

  nup=nu0*(1.0+zp)/(1.0+z);

  dJdz=SPEEDOFLIGHT_CGS/4.0/PI;
  dJdz *= (1.0+z)*(1.0+z)*a->radioEmission(nup,zp);
  dJdz /= c->hubbleZ(zp)*UNH*c->getH();
  
  return dJdz;
}

double getDJradioDz(double zp)
{
  return setDJradioDz(0.0,0.0,zp,NULL,NULL,0);
}


//////////////////////////////////////////////////////////////////////////
//  Radio Loud AGN model
/////////////////////////////////////////////////////////////////////////
// Haiman, Quataert, Bower (2004)  ApJ 612, 698

// Calculate central black hole mass within halo
// based upon M-sigma relation
//
// Units: mh   M_sol
//        mbh  M_sol
//
double Astrophysics::blackHoleMass(double mh, double z)
{
  double mbh;

  mbh=1.0e6*pow(mh/1.5e12,5.0/3.0);
  mbh*=pow(c->getH()/0.71,5.0/3.0);
  mbh*=pow(1.0+z,2.5);

  return mbh;
}

// Calculate eddington luminosity
//
// Units: mb  Msol
//        lum Erg s^{-1}
//
double Astrophysics::luminosityEddington(double mb)
{
  double lum;

  lum=1.0e46*(mb/1.0e8);
  return lum;
}

// Calculate radio loudness
// Units:  lumR  Erg s^{-1} Hz^{-1}
//
double Astrophysics::luminosityRadioAGN(double nu,double mh, double z, double R)
{
  double lumR,lumE;
  double mb;
  double alphaR(0.0);
  double xrayfrac(4.0e-16);  //this value is fudged to match HQB: CHECK IT!!!

  mb=blackHoleMass(mh,z);
  lumE=luminosityEddington(mb);

  lumR=lumE*xrayfrac;   // correct for fraction of energy in Xray
  lumR*=pow(10.0,R);  //fraction in radio flux
  lumR*=pow(nu/2.8e9,-alphaR);

  return lumR;
}

// Calculate fractional number of radio loud sources
//
double Astrophysics::numberRadioAGN(double R)
{
  double nR;
  double fracLoud(0.1);
  double sig2(4.0/PI);
  double Rbar(2.8);

  nR=0.5*fracLoud*exp(-pow(R-Rbar,2.0)/sig2);

  return nR;
}

//Calculate duty cycle of radio loud quasars
//
double Astrophysics::dutycycleRadioAGN(double z)
{
  double tH;
  double tQ;

  tQ=2.0e7*YEAR;
  tH=c->cosmicTime(z)/UNH/c->getH();
 
  return tQ/tH;
}


//Get observed flux from source luminosity
//
//  Units:  mh  msol
//          fluxR  Erg cm^{-2} s^{-1} Hz^{-1}
//
double Astrophysics::fluxRadioAGN(double nu, double mh, double z, double R)
{
  double lumR, dL;
  double zz;
  double fluxR;

  zz=z+1.0;

  lumR=luminosityRadioAGN(nu*zz,mh,z,R);
  dL=c->lumDistance(z)*SPEEDOFLIGHT_CGS/UNH/c->getH();
  fluxR=lumR/(4.0*PI*dL*dL);

  return fluxR;
}

// now get number statistics for all of this

// Find R required for radio flux to exceed flux limit Flim
// given halo mass mh at redshift z
//
double Astrophysics::getRfromFluxLimit(double nu, double mh, double z, double Flim)
{
  double R;
  double rat;

  rat=Flim/fluxRadioAGN(nu,mh,z,0.0);  //rat=10^R
  R=log(rat)/log(10.0);
  return R;
}

// Find N(R>Rlim)
//
double Astrophysics::fracBrighterR(double nu, double mh, double z, double Flim)
{
  double Rlim;
  double fracR;
  double tol(1.0e-4);
  double Rmax(5.0e2);  // need to handle this integral better
  double mb;
  double mbcut(0.0);   // use to restrict AGN to large black holes

  mb=blackHoleMass(mh,z);
  if(mb<mbcut) return 0.0;

  Rlim=getRfromFluxLimit(nu,mh,z,Flim);

  setDummyNumberRadioAGN(0.0,this,1);

  fracR=qromb(dummyNumberRadioAGN,Rlim,Rmax,tol);   // integrate from Rlim to inf
  return fracR;
}

double dummyNumberRadioAGN(double R)
{
  return setDummyNumberRadioAGN(R,NULL,0);
}

double setDummyNumberRadioAGN(double R, Astrophysics *a1, int iflag)
{
  static Astrophysics *a;

  if(iflag==1){
    a=a1;
    return 0.0;
  }

  return a->numberRadioAGN(R);
}


// Calculate number density of objects with F>Flim
// Uses Press-Schecter mass function
//
//
double Astrophysics::numberDensRadioAGNBrighterF(double nu, double z, double Flim)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double number;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setRadioNumberIntegrand(0.0,ans,ans,nu,Flim,z,c,this,1);
  mMin=coolMass(c,z);  //minimum halo mass

  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];

    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
    	   radioNumberIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6)
      break;
    mMin = mMax;
  }
  number = ans[1];
  //number /= c->fColl(z);  //normalise halo mass function appropriately
  free_dvector(ans,1,1);

  return number;

}

void setRadioNumberIntegrand(double m, double y[], double deriv[],double nu1, 
			     double Flim1, double z1, Cosmology *c1, Astrophysics *a1,int flag)
{
  static double Flim;
  static double nu;
  static double z;
  static Cosmology *c;
  static Astrophysics *a;
  double weight;
  double mb;
  double mbcut(1.0e7);
  double lowfrac(1.0e-7);

  if (flag == 1) {
    nu=nu1;
    Flim=Flim1;
    z = z1;
    c=c1;
    a=a1;
    return;
  }


  mb=a->blackHoleMass(m,z);
  weight=a->fracBrighterR(nu,m,z,Flim);
  if(mb<mbcut)  weight*=lowfrac;
  //cout<<weight<<endl;
  deriv[1] = nmST(m,z,c)/m*weight;
  //  deriv[1] = c->dndlMJenkins(z,m)/m*weight;
  
}

// Dummy function for integrating the bubble filling Q
void radioNumberIntegrand(double m, double y[], double deriv[])
{
  setRadioNumberIntegrand(m,y,deriv,0.0,0.0,0.0,NULL,NULL,0);
}

// Calculate number density per unit redshift and solid angle of radio
// loud AGN with flux greater than Flim.
//  dN/dzdOmega
//
double Astrophysics::numberRadioAGNBrighterF(double nu, double z, double Flim)
{
  double nR;
  double NRadio;
  double volume;
  double cH,r;

  nR=numberDensRadioAGNBrighterF(nu,z,Flim);  //comoving number density
  cH=SPEEDOFLIGHT_CGS/c->getH()/H0_CMSMPC;
  r=cH*c->coordDistance(z);
  volume=cH/c->hubbleZ(z)*r*r;    //comoving volume element

  NRadio=dutycycleRadioAGN(z)*nR*volume;

  return NRadio;
}

//Number counts of radio loud AGN
double Astrophysics::numberCountsRadioAGN(double Flim, double zMin, double zMax, double nu)
{
  double tol(1.0e-4);
  double counts;
 
  if(zMax<0.0) zMax=ZSFR;
  
  setNumberCountsAGNInt(0.0,nu,Flim,this,1);
  counts=qromb(numberCountsAGNInt,zMin,zMax,tol);
    
  return counts;
}

double setNumberCountsAGNInt(double z, double nu1, double Flim1, Astrophysics *a1, int iflag)
{
  static double nu;
  static double Flim;
  static Astrophysics *a;

  if(iflag==1){
    nu=nu1;
    Flim=Flim1;
    a=a1;
    return 0.0;
  }
  
  return a->numberRadioAGNBrighterF(nu,z,Flim);
}

double numberCountsAGNInt(double z)
{
  return setNumberCountsAGNInt(z,0.0,0.0,NULL,0);
}


// Calculate comoving emissivity of radio loud AGN
// Uses Press-Schecter mass function
//
// Units:  emiss  Erg s^-1 Hz^-1 Mpc^-3
//
double Astrophysics::emissivityRadioAGN(double nu, double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double emiss;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setRadioEmissIntegrand(0.0,ans,ans,nu,z,c,this,1);
  mMin=coolMass(c,z);  //minimum halo mass

  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];

    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
    	   radioEmissIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6)
      break;
    mMin = mMax;
  }
  emiss = ans[1];
  //number /= c->fColl(z);  //normalise halo mass function appropriately
  free_dvector(ans,1,1);

  return emiss;
}

void setRadioEmissIntegrand(double m, double y[], double deriv[],double nu1, 
			     double z1, Cosmology *c1, Astrophysics *a1,int flag)
{
  static double nu;
  static double z;
  static Cosmology *c;
  static Astrophysics *a;
  double weight;
  double R(2.8);
  double fracLoud(0.1);
  double mb;
  double mbcut(1.0e7);
  double lowfrac(1.0e-7);
  double lowR(1.6);

  if (flag == 1) {
    nu=nu1;
    z = z1;
    c=c1;
    a=a1;
    return;
  }

  mb=a->blackHoleMass(m,z);
  if(mb<mbcut){
    //    R=lowR;
    weight=a->luminosityRadioAGN(nu,m,z,R)*fracLoud*lowfrac;
  } else{
    weight=a->luminosityRadioAGN(nu,m,z,R)*fracLoud;
  }
  weight*=a->dutycycleRadioAGN(z);

  deriv[1] = nmST(m,z,c)/m*weight;
  
}

// Dummy function for integrating the bubble filling Q
void radioEmissIntegrand(double m, double y[], double deriv[])
{
  setRadioEmissIntegrand(m,y,deriv,0.0,0.0,NULL,NULL,0);
}

double Astrophysics::differentialNumberCountsAGN(double S, double nu, double zobs)
{
  double dNdS;
  double zmin,zmax;
  double tol(1.0e-4);
  
  zmin=zobs;
  zmax=ZSFR;

  setDiffNumberCountAGNInt(zobs,nu,S,this,c,1);

  dNdS=qromb(diffNumberCountAGNInt,zmin,zmax,tol);

  return dNdS;
}

double setDiffNumberCountAGNInt(double z, double nu1, double S1, Astrophysics *a1, Cosmology *c1, int iflag)
{
  static double nu;
  static Astrophysics *a;
  static Cosmology *c;
  static double S;
  static double zobs;
  double mh;
  double dmds;
  double dndsdv,volume;
  double mCool;
  double fracloud(0.1);
  double R(2.8);

  if(iflag==1){
    nu=nu1;
    a=a1;
    c=c1;
    S=S1;
    zobs=0.0;
    return 0.0;
  }


  dndsdv=a->dndsdvRadioAGN(S,nu,z);
  volume=c->volumeComoving(z);

  return dndsdv*volume;
}

double diffNumberCountAGNInt(double z)
{
  return setDiffNumberCountAGNInt(z,0.0,0.0,NULL,NULL,0);
}

// Calculate differential number density of radio loud AGN
// Uses Press-Schecter mass function
//
// Units:  dndsdv (Erg s^-1 Hz^-1)^-1 Mpc^-3
//
double Astrophysics::dndsdvRadioAGN(double S, double nu, double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double dndsdv;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setRadioAGNdNdSIntegrand(0.0,ans,ans,S,nu,z,c,this,1);
  mMin=coolMass(c,z);  //minimum halo mass

  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];

    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
    	   radioAGNdNdSIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6)
      break;
    mMin = mMax;
  }
  dndsdv = ans[1];
  //number /= c->fColl(z);  //normalise halo mass function appropriately
  free_dvector(ans,1,1);

  return dndsdv;
}

void setRadioAGNdNdSIntegrand(double m, double y[], double deriv[],double S1, double nu1, double z1, Cosmology *c1, Astrophysics *a1,int flag)
{
  static double nu;
  static double z;
  static Cosmology *c;
  static Astrophysics *a;
  static double S;
  double weight,drds;
  double R,ratio;
  double NR;
  //  double mbcut(1.0e7);
  //double mb;
  //double small_frac(1.0e-7);

  if (flag == 1) {
    S=S1;
    nu=nu1;
    z = z1;
    c=c1;
    a=a1;
    return;
  }

  //  mb=a->blackHoleMass(m,z);

  ratio=S/a->fluxRadioAGN(nu,m,z,0.0);
  R=log(ratio)/log(10.0);

  drds=1.0/a->fluxRadioAGN(nu,m,z,R)/log(10.0);
  NR=a->numberRadioAGN(R);

  weight=a->dutycycleRadioAGN(z)*NR*drds*a->massCutRadioAGN(m,z);
  //  cout<<m<<"\t"<<z<<"\t"<<R<<"\t"<<NR<<"\t"<<dmds<<"\t"<<weight<<endl;
  //  if(mb<mbcut) weight*=small_frac;
  deriv[1] = nmST(m,z,c)/m*weight;
  
}

// Dummy function for integrating the bubble filling Q
void radioAGNdNdSIntegrand(double m, double y[], double deriv[])
{
  setRadioAGNdNdSIntegrand(m,y,deriv,0.0,0.0,0.0,NULL,NULL,0);
}

double Astrophysics::massCutRadioAGN(double mh, double z)
{
  double lowfrac(1.0e-7);
  double mbcut(1.0e7);
  double mb;
  double weight;

  mb=blackHoleMass(mh,z);
  if(mb<mbcut) return lowfrac;
  return 1.0;
}


////////////////////////////////////////////////////////////////////////
// Number count code
////////////////////////////////////////////////////////////////////////

double Astrophysics::fluxMoments(double sMin, double sMax, double zMin, double zMax, double nu, int order)
{
  double tol(1.0e-4);
  double moment;
 
  if(zMax<0.0) zMax=ZSFR;
  
  setFluxMomentsInt(0.0,nu,sMin,sMax,order,this,c,1);
  moment=qromb(fluxMomentsInt,zMin,zMax,tol);
    
  return moment;
}

double setFluxMomentsInt(double z, double nu1, double sMin1, double sMax1, int order1, Astrophysics *a1, Cosmology *c1,int iflag)
{
  static double nu;
  static double sMin;
  static double sMax;
  static int order;
  static Astrophysics *a;
  static Cosmology *c;
  double volume;

  if(iflag==1){
    nu=nu1;
    sMin=sMin1;
    sMax=sMax1;
    order=order1;
    a=a1;
    c=c1;
    return 0.0;
  }

  volume=c->volumeComoving(z);
  return a->fluxMomentsDmDz(z,nu,sMin,sMax,order)*volume;
}


double fluxMomentsInt(double z)
{
  return setFluxMomentsInt(z,0.0,0.0,0.0,0,NULL,NULL,0);
}


// Calculate d<S^n>/dz
//
//
double Astrophysics::fluxMomentsDmDz(double z, double nu, double sMin1, double sMax1, int order)
{
  double sMin,sMax;
  double mMin, mMax, mEnd,mCool;
  double mINFTY(1.0e90);
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double number;
  int b(0);

  sMin=sMin1;
  sMax=sMax1;

  mMin=massFlux(sMin,z,nu);
  mEnd=massFlux(sMax,z,nu);
  mCool=coolMass(c,z);

  if(sMax1<0.0){  // max(S_cut, S_cool) to \infty
    mEnd=mINFTY;
    mMin=fmax(mMin,mCool);
  }else{         // S_cool to S_cut
    mMin=mCool;
  }

  //perform integration
  ans = dvector(1,1);
  ans[1] = 0.0;
  setFluxMomentIntegrand(0.0,ans,ans,nu,z,order,c,this,1);

  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    if(mMax>mEnd){
      mMax=mEnd;
      b=1;
    }

    oldans = ans[1];

    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
    	   fluxMomentIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4 || b==1)
      break;
    mMin = mMax;
  }
  number = ans[1];
  //  number /= CRITDENMSOLMPC*c->getOm0hh(); //for testing fcoll agreement

  free_dvector(ans,1,1);

  return number;
}


void setFluxMomentIntegrand(double m, double y[], double deriv[],double nu1, 
			     double z1, int order,Cosmology *c1, Astrophysics *a1,int flag)
{
  static double nu;
  static double z;
  static double nOrder;
  static Cosmology *c;
  static Astrophysics *a;
  double weight;
  double S,dndsdv;
  double tSB(1.0e7);  //starburst duration in years
  double htSB;
  static double zSB;

  if (flag == 1) {
    nu=nu1;
    z = z1;
    c=c1;
    a=a1;
    nOrder=(double)(order);

    htSB=c->cosmicTime(z)+tSB*YEAR*UNH*c->getH();
    zSB=c->zFromT(htSB,100.0);
    return;
  }

  S=a->fluxMass(m,z,nu);

  dndsdv=a->haloFluxDensity(zSB,z,m,nu);

  weight=dndsdv*a->fluxMass(m,z,nu)/m;
  weight*=pow(S,nOrder);
  //  weight=nmDZ(m,z,c); //for testing fcoll agreement

  deriv[1] = weight;
  
}

// Dummy function for integrating the bubble filling Q
void fluxMomentIntegrand(double m, double y[], double deriv[])
{
  setFluxMomentIntegrand(m,y,deriv,0.0,0.0,0,NULL,NULL,0);
}

// for the moment just use Free-free flux
// Units:  flux: ergs s^-1 Hz^-1
double Astrophysics::fluxMass(double mass, double z, double nu)
{
  double lum, flux,lD,cH;

  cH=SPEEDOFLIGHT_CGS/c->getH()/UNH;
  lD=c->lumDistance(z)*cH;

  //figures direct from Oh1999 paper
  //  lum=2.0*(mass/1.0e9)*(FSTAR/0.17);
  //lum*=1.2e27; //ff

  //  cout<<lum<<"\t"<<lumRadioFreeFree(0.0,haloStarFormationRate(mass))/lum<<endl;

  //  lum=lumRadioFreeFree(nu*(1.0+z),haloStarFormationRate(mass));
  lum=lumRadioGalaxy(nu*(1.0+z),haloStarFormationRate(mass));

  //  lum*=1.4e41;  //Halpha
  flux=lum/4.0/PI/lD/lD;
  flux*=1.0+z;   //flux per unit frequency not total flux

  return flux;
}

//assumes flux scales linearly with mass
double Astrophysics::massFlux(double flux, double z, double nu)
{
  double mass;

  mass=flux/fluxMass(1.0e9,z,nu)*1.0e9;

  return mass;
}

// Specific Luminosity of a radio galaxy
// Units: sE : ergs s^-1 Hz^-1
//        sfr : Msol Year^-1
double Astrophysics::lumRadioGalaxy(double nu, double sfr)
{
  double sE;
  double emission0(1.5e28);
  double spec_alpha(-0.8);

  sE=0.06/1.44*pow(nu/1.49e9,-0.5);
  sE+=1.38/1.44*pow(nu/1.49e9,spec_alpha);
  sE*=emission0;
  sE*=sfr;

  return sE;
}

//specific luminosity of free-free ionized halo
// Units:   ergs s^-1 Hz^-1
//          sfr Msol Year^-1
double Astrophysics::lumRadioFreeFree(double nu, double sfr)
{
  double lum(1.2e27);
  //should have some dependence on NION in here - need pivot value
  lum*=sfr*(1.0-FESC);
  return lum;
}

//star formation rate of a halo
// assumes the Haiman and Loeb 1998 prescription
double Astrophysics::haloStarFormationRate(double mh)
{
  double mg;
  double sfr;

  mg=c->getOmegab()*mh/c->getOmega0();
  sfr=FSTAR*mg/STARBURSTAGE;
  return sfr;
}

// Calculate number of halos in flux interval given
// redshift z and starburst duration corresponds to redshift zSB
//
double Astrophysics::haloFluxDensity(double z, double zSB, double M, double nu)
{
  double dndsdv;
  double tol(1.0e-4);

  setHaloFluxDensityInt(z,M,nu,this,c,1);

  dndsdv=qromb(haloFluxDensityInt,z,zSB,tol);

  return dndsdv;
}

double setHaloFluxDensityInt(double z, double mh1, double nu1, Astrophysics *a1,Cosmology *c1, int iflag)
{
  static double mh;
  static Cosmology *c;
  static Astrophysics *a;
  static double nu;
  double weight;

  if(iflag==1){
    nu=nu1;
    mh=mh1;
    c=c1;
    a=a1;
    return 0.0;
  }

  //slight redshift dependence in S(M)
  weight=mh/a->fluxMass(mh,z,nu);   //dMdS - linear mass dependence assumed

  return nmDZ(mh,z,c)/mh*weight;
}

double haloFluxDensityInt(double z)
{
  return setHaloFluxDensityInt(z,0.0,0.0,NULL,NULL,0);
}

double Astrophysics::differentialNumberCountsRG(double S, double nu, double zobs)
{
  double dNdS;
  double zmin,zmax;
  double tol(1.0e-4);
  
  zmin=zobs;
  zmax=ZSFR;

  setDiffNumberCountRGInt(zobs,nu,S,this,c,1);

  dNdS=qromb(diffNumberCountRGInt,zmin,zmax,tol);

  return dNdS;
}

double setDiffNumberCountRGInt(double z, double nu1, double S1, Astrophysics *a1, Cosmology *c1, int iflag)
{
  static double nu;
  static Astrophysics *a;
  static Cosmology *c;
  static double S;
  static double zobs;
  double m;
  double htSB,zSB;
  double dndsdv,volume;
  double tSB(1.0e7);
  double mCool;

  if(iflag==1){
    nu=nu1;
    a=a1;
    c=c1;
    S=S1;
    zobs=0.0;
    return 0.0;
  }

  m=a->massFlux(S,z,nu);
  //  mCool=coolMass(c,z);
  //if(m<mCool) return 0.0;

  htSB=c->cosmicTime(z)+tSB*YEAR*UNH*c->getH();
  zSB=c->zFromT(htSB,100.0);
  dndsdv=a->haloFluxDensity(zSB,z,m,nu);
  volume=c->volumeComoving(z);

  return dndsdv*volume;
}

double diffNumberCountRGInt(double z)
{
  return setDiffNumberCountRGInt(z,0.0,0.0,NULL,NULL,0);
}
