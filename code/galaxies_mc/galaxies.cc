

/*  galaxies.cc
 *
 *  Stuff for calculating emissivity from high redshift galaxies
 *
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include "astrophysics.h"
#include "dcosmology.h"
#include "dnumrecipes.h"
#include "galaxies.h"
#include "haloDensity.h"
#include "reionization.h"
#include <sstream>
#include <gsl/gsl_sf.h>


using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

Galaxies::Galaxies()
{

}

Galaxies::~Galaxies()
{
  // cout <<"Reionization destructor has been called." <<endl;

}

////////////////////////////////////////////////////////////////////
// Initiation routines
////////////////////////////////////////////////////////////////////
void Galaxies::initGalaxies(Cosmology *c1, Astrophysics *a1)
{
  c=c1;
  a=a1;
}

void Galaxies::setCosmology(Cosmology *c1)
{
   c=c1;
}

////////////////////////////////////////////////////////////////////
// Galaxy luminosity function
////////////////////////////////////////////////////////////////////
//Calculate galaxy luminosity function using the Schecter form 
//- Schechter (1976)
// three free parameters for the fit:
//   phi0 : Mpc^{-3}
//   alpha 
//   L0   : ergs s^{-1}  or erg s^{-1} Hz^{-1}
// units: phi: no. galaxies per unit comoving volume per unit luminosity
double Galaxies::luminosityFunction(double L, double alpha, double phi0, double L0)
{
  double phi;
  phi=phi0/L0;
  phi*=pow(L/L0,alpha);
  phi*=exp(-L/L0);
  return phi;
}

//number density of galaxies from integrating the Schecter function above some 
//minimum luminosity L0
double Galaxies::numberDensityLF(double Lmin, double alpha, double phi0, double L0)
{
   //cout<<phi0<<"\t"<<1.0+alpha<<"\t"<<Lmin<<"\t"<<L0<<endl;
  return phi0*gsl_sf_gamma_inc(1.0+alpha,Lmin/L0);
}

//luminosity density from integrating the Schecter function above some 
//minimum luminosity L0
double Galaxies::lumDensityLF(double Lmin, double alpha, double phi0, double L0)
{
  return phi0*L0*gsl_sf_gamma_inc(2.0+alpha,Lmin/L0);
}

//Calculate galaxy luminosity function using the Schecter form 
// in terms of absolute magnitude M
// three free parameters for the fit:
//   phi0 : Mpc^{-3}
//   alpha 
//   M0   :
// units: phi: no. galaxies per unit comoving volume per unit luminosity
//
// Note: dn/dL very different from dn/dM or dn/dm due to units
double Galaxies::luminosityFunctionM(double M, double alpha, double phi0, double M0)
{
  double phi;
  phi=phi0*log(10.0)/2.5;
  phi*=exp(-0.4*(M-M0)*(alpha+1.0)*log(10.0));
  phi*=exp(-exp(-0.4*(M-M0)*log(10.0)));

  return phi;
}

//number density of galaxies from integrating the Schecter function above some 
//minimum luminosity L0
double Galaxies::numberDensityLF_M(double Mmin, double alpha, double phi0, double M0)
{
   double Lmin, L0,ngal;
   double z(1.0); //cancels out
   Lmin=lumFromAbsMag(Mmin,z);
   L0=lumFromAbsMag(M0,z);
   ngal= phi0*gsl_sf_gamma_inc(1.0+alpha,Lmin/L0);
   //cout<<Lmin<<"\t"<<L0<<"\t"<<ngal<<endl;
   return ngal;
}

//use fit from Bouwyens+ 2008 ApJ 686, 230
//
//Units: phi Mpc^-3
//would be nice to output errors too
vector<double> Galaxies::getLFParam(double z)
{
   double M0, sig_M0, M1, sig_M1;
   double phi0, sig_phi0, phi1, sig_phi1;
   double alpha0, sig_alpha0, alpha1, sig_alpha1;
   double z0(3.8);
   double M, phi, alpha;
   vector<double> param;

   M0=-21.02;
   sig_M0=0.09;
   M1=0.36;
   sig_M1=0.08;
   
   phi0=1.16;
   sig_phi0=0.20;
   phi1=0.024;
   sig_phi1=0.065;
   
   alpha0=-1.74;
   sig_alpha0=0.05;
   alpha1=0.04;
   sig_alpha1=0.05;

   M=M0+M1*(z-z0);
   phi=phi0*pow(10.0,phi1*(z-z0))*1.0e-3;
   alpha=alpha0+alpha1*(z-z0);

   param.push_back(M);
   param.push_back(phi);
   param.push_back(alpha);   
   return param;
}



////////////////////////////////////////////////////////////////////
// AB magnitude conversions
////////////////////////////////////////////////////////////////////

// Calculates the apparent AB magnitude given the luminosity
//Units: L  erg s^-1 Hz^-1
//       flux erg s^-2 cm^-2 Hz^-1
double Galaxies::magFromL(double L, double z)
{
  double flux;
  double dL, cH;
  double m;
  cH=SPEEDOFLIGHT_CGS/(UNH*c->getH());
  dL=c->lumDistance(z)*cH;

  //calculate flux assuming luminosity in a freq. interval
  flux=L*(1.0+z)/(4.0*PI*dL*dL);
  m= -2.5*log10(flux)-48.60;

  return m;
}

// Calculates the luminosity given the apparent AB magnitude
//Units: L  erg s^-1 Hz^-1
double Galaxies::lumFromMag(double m, double z)
{
  double flux;
  double dL, cH;
  double L;
  cH=SPEEDOFLIGHT_CGS/(UNH*c->getH());
  dL=c->lumDistance(z)*cH;
  L=exp(-0.4*(m+48.60)*log(10.0));
  L*=4.0*PI*dL*dL/(1.0+z);

  return L;
}

// Calculates the absolute AB magnitude given the luminosity
//Units: L  erg s^-1 Hz^-1
double Galaxies::absMagFromL(double L, double z)
{
  double M;
  double dL, cH;
  cH=SPEEDOFLIGHT_CGS/(UNH*c->getH());
  dL=c->lumDistance(z)*cH*1.0e6/MPC;
  M=magFromL(L,z)+2.5*log10(1+z)-5.0*log10(dL)+5.0;
  return M;
}

// Calculates the absolute AB magnitude given the relative mag
//Units: L  erg s^-1 Hz^-1
//Note: expects dL in Parsecs
double Galaxies::absMagFromM(double m, double z)
{
  double M;
  double dL, cH;
  cH=SPEEDOFLIGHT_CGS/(UNH*c->getH());
  dL=c->lumDistance(z)*cH*1.0e6/MPC;
  M=m+2.5*log10(1+z)-5.0*log10(dL)+5.0;
  return M;
}

// Calculates the luminosity given the absolute AB magnitude
//Units: L  erg s^-1 Hz^-1
//
// Note: expects dL in Parsecs
double Galaxies::lumFromAbsMag(double M, double z)
{
  double m, L;
  double dL, cH;

  //convert M to apparent magnitude
  cH=SPEEDOFLIGHT_CGS/(UNH*c->getH());
  dL=c->lumDistance(z)*cH*1.0e6/MPC;
  m=M-2.5*log10(1+z)+5.0*log10(dL)-5.0;
//NOTE think I may be confused over luminosity and 
// luminosity/Hz in this expressions

  //use conversion from apparent mag to L
  L=lumFromMag(m,z);

  cout<<"L1="<<L<<endl;
  //simpler conversion
  //M= -2.5*log10(L)+5.0+2.5*log10(4.0*PI)-48.60;
  
  L= -M-48.60+5.0+2.5*log10(4.0*PI);
  L=pow(10.0,L/2.5);
  cout<<"L2="<<L<<endl;

  return L;
}

////////////////////////////////////////////////////////////////////
// Galaxy contribution to reionization
////////////////////////////////////////////////////////////////////
//Calculate the ionizing photon emissivity at the Lyman limit
// in usints of 10^25 erg s^-1 Mpc^-1
//
// Lmin = lower bound on luminosty
// nuObs is frequency at which luminosity function is measured
// phi0, L0, alpha describe the galaxy luminosity function of the sample
//
// This is epsilon_25 defined in Bolton & Haehnelt 2007
double Galaxies::emissivityGal(double Lmin, double nuObs, double z, double phi0, double L0, double alpha)
{
  double epsilon;

  epsilon=lumDensityLF(Lmin,alpha,phi0,L0);
  epsilon/=nuObs;
  epsilon/=6.0;  //break between lambda=1500A and lyman edge
  epsilon/=1.0e25;
  return epsilon;
}

// rate of production of ionizing photons per unit volume
// ndot s^-1 Mpc^-3
//
double Galaxies::ndotGal(double e25, double alphaS, double fesc)
{
  double ndot;
  double ndot0(5.03063e49*e25);
  ndot=ndot0*(fesc/0.1)/(alphaS/3.0);

  return ndot;
}

// rate of production of ionizing photons per unit volume
// ndot s^-1 Mpc^-3
// mfp  Mpc
//
//
double Galaxies::gammaGal(double e25, double alphaS, double fesc, double mfp, double z)
{
  double gamma;
  double gamma0(0.0326199);
  gamma=gamma0*e25*(fesc/0.1)/((alphaS+3.0)/6.0)*(mfp/40.0);
  gamma*=pow((1.0+z)/7.0,2.0);

  return gamma;
}

// Rate of recombinations needed to just maintain constant ionized fraction
// ndot s^-1 Mpc^-3
//
double Galaxies::recombLevel(double z, double xi, double C)
{
   double recomb;
   double nh;
   double alphaB;
   double tg(1.0e4);

   nh=c->nh(0.0);
   alphaB=a->recombH(tg);
   recomb=xi*C*nh*nh*alphaB*pow(1.0+z,3.0);
   recomb*=MPC*MPC*MPC;  //convert from cm^-3 to Mpc^-3

  return recomb;
}

////////////////////////////////////////////////////////////////////
// Quasar contribution to reionization
////////////////////////////////////////////////////////////////////
//Calculate the ionizing photon emissivity at the Lyman limit
// in usints of 10^24 erg s^-1 Mpc^-1
//
// Lmin = lower bound on luminosty
// nuObs is frequency at which luminosity function is measured
// phi0, L0, alpha describe the galaxy luminosity function of the sample
//
// This is epsilon_24 defined in Bolton & Haehnelt 2007
double Galaxies::emissivityQSO(double Lmin, double nuObs, double z, double phi0, double L0, double beta1, double beta2)
{
  double epsilon;
  double lamObs;

  epsilon=lumDensityLF_QSO(Lmin,phi0,L0,beta1,beta2);
  epsilon/=nuObs;

  //connect to 912 angstrom with a brokwn power law
  lamObs=SPEEDOFLIGHT_CGS/nuObs;
  epsilon*=pow(1050.0e-8/lamObs,-0.5);
  epsilon*=pow(912.0/1050.0,-1.5);

  epsilon/=1.0e24;
  return epsilon;
}

//Calculate quasar luminosity function using a double power law form
// d\Phi/dL   (7) in Hopkins et al (2007)
// four free parameters for the fit:
//   phi0 : Mpc^{-3}
//   beta1, beta2 
//   L0   : ergs s^{-1}  or erg s^{-1} Hz^{-1}
// units: phi: no. galaxies per unit comoving volume per unit luminosity
double Galaxies::luminosityFunctionQSO(double L, double phi0, double L0, double beta1, double beta2)
{
  double phi;
  phi=phi0/L0;
  phi/=pow(L/L0,-beta1)+pow(L/L0,-beta2);
  return phi;
}

//Calculate quasar luminosity function using a double power law form
// d\Phi/dlog L   (6) in Hopkins et al (2007)
// four free parameters for the fit:
//   phistar : Mpc^{-3}
//   beta1, beta2 
//   L0   : ergs s^{-1}  or erg s^{-1} Hz^{-1}
// units: phi: no. galaxies per unit comoving volume per unit luminosity
double Galaxies::luminosityFunctionQSO_dNdLogL(double L, double phistar, double L0, double gamma1, double gamma2)
{
  double phi;
  phi=phistar;
  phi/=pow(L/L0,gamma1)+pow(L/L0,gamma2);
  return phi;
}

//Calculate quasar luminosity function using a double power law form
// d\Phi/dlog L   (7) in Hopkins et al (2007)
// four free parameters for the fit:
//   phistarp = : Mpc^{-3}
//   beta1, beta2 
//   L0   : ergs s^{-1}  or erg s^{-1} Hz^{-1}
// units: phi: no. galaxies per unit comoving volume per unit luminosity
// phistarp = phistar/ ln10
// alpha= -(gamma1+1) ; beta= -(gamma2+1)
double Galaxies::luminosityFunctionQSO_dNdL(double L, double phistarp, double L0, double alpha, double beta)
{
  double phi;
  phi=phistarp/L0;
  phi/=pow(L/L0,-alpha)+pow(L/L0,-beta);
  return phi;
}

//Calculate quasar luminosity function using a double power law form
// d\Phi/dM   (8) in Hopkins et al (2007)
// four free parameters for the fit:
//   phistarpp : Mpc^{-3}
//   alpha, beta 
//   M0   : absolute magnitude
// units: phi: no. galaxies per unit comoving volume per unit luminosity
// phistarpp=0.4*phistar
double Galaxies::luminosityFunctionQSO_dNdM(double M, double phistarpp, double M0, double alpha, double beta)
{
  double phi;
  phi=pow(10.0,0.4*(alpha+1.0)*(M-M0));
  phi+=pow(10.0,0.4*(beta+1.0)*(M-M0));
  phi=phistarpp/phi;
  return phi;
}


//luminosity density from integrating the Schecter function above some 
//minimum luminosity Lmin
//
// this is an integration done as \int dL dn/dL
double Galaxies::lumDensityLF_QSO(double Lmin, double phi0, double L0, double beta1, double beta2)
{
   double tol(1.0e-4);
   double xmin, xmax;

   xmin=Lmin/L0;
   xmax=10.0;   //stand in for infinity - should drop rapidly above L0
   setDummyLFQSO(0.0,beta1,beta2,this,1);
   return phi0*L0*qromb(dummyLFQSO,xmin,xmax,tol);
}

double dummyLFQSO(double x)
{
   return setDummyLFQSO(x,0.0,0.0,NULL,0);
}

double setDummyLFQSO(double L, double beta1in, double beta2in, Galaxies *Gal1, int iflag)
{
   static double beta1, beta2;
   static Galaxies *Gal;
   double phi;

   if(iflag==1){
      beta1=beta1in;
      beta2=beta2in;
      Gal=Gal1;
      return 0.0;
   }

   phi=Gal->luminosityFunctionQSO(L,1.0,1.0,beta1,beta2);
   return phi;
}

// rate of production of ionizing photons per unit volume
// ndot s^-1 Mpc^-3
//
double Galaxies::ndotQSO(double e24, double alphaS)
{
  double ndot;
  double ndot0(5.03063e49*e24);
  ndot=ndot0/(alphaS/3.0);

  return ndot;
}

// rate of production of ionizing photons per unit volume
// ndot s^-1 Mpc^-3
// mfp  Mpc
//
//
double Galaxies::gammaQSO(double e24, double alphaS, double mfp, double z)
{
  double gamma;
  double gamma0(0.04);
  gamma=gamma0*e24/((alphaS+3.0)/4.5)*(mfp/40.0);
  gamma*=pow((1.0+z)/7.0,2.0);

  return gamma;
}

///////////////////////////////////////////////////////
// Hopkins et al LDDE QSO luminosity function
//////////////////////////////////////////////////////
//Return QSO luminosity function
//using fit from Hopkins, Gordon, Hernquist (2007) - ApJ 654, 731
// for Luminosity-dependent Density Evolution (LDDE)
//from Eq. (11)-(10) + Table 4
//
// Expects Luminosty L [erg s^{-1}]
//
// Return: [log10(phi0/Mpc^{-3}),Log10(L0/Lsol),gamma1,gamma2]
// where Lsol=3.9e33 erg s^{-1}
//
//Units: phi Mpc^-3
//
double Galaxies::luminosityfunctionQSO_LDDE(double L, double z)
{
   double phi,phistar, Lstar;
   double lphi0, lL0, gamma1,gamma2;
   vector<double> param;
   double lLstar0,lLc, zc0,alpha;
   double p1_46, p2_46;
   double beta1, beta2;
   double rhoBH0;
   double ed;
   double p1,p2,Lc,zc;

   lphi0= -6.2;
   lLstar0= 45.99;
   gamma1= 0.933;
   gamma2= 2.20;
   lLc= 46.72;
   zc0=1.852;
   alpha=0.274;
   p1_46=5.95;
   p2_46= -1.65;
   beta1=0.29;
   beta2= -0.62;
   rhoBH0=5.09;


   //powerlaw indices
   p1=p1_46+beta1*log10(L/1.0e46);
   p2=p2_46+beta2*log10(L/1.0e46);

   //threshold redshift
   if(L<=Lc){
      zc=zc0*pow(L/Lc,alpha);
   }else{
      zc=zc0;
   }

   //fitting function e_d(L,z), which provides multiplicative 
   //envelope for the normal Luminosity function
   if(z<=zc){
      ed=pow(1.0+z,p1);
   }else{
      ed=pow(1.0+zc,p1)*pow((1.0+z)/(1.0+zc),p2);
   }

   //get parameters in 
   phistar=pow(10.0,lphi0);
   Lstar=pow(10.0,lL0);

   //prepare phi(L)
   phi=phistar;
   phi/=pow(L/Lstar,gamma1)+pow(L/Lstar,gamma2);
   phi*=ed;

   return phi;
}

/////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////
//number density from integrating the Schecter function above some 
//minimum luminosity Lmin in Lsol
// at a redshift z
// using the model iflag = 0: full  1: PLE
//
// this is an integration done as \int dL dn/dL
double Galaxies::numberDensityQSO_HOP(double z, double Lmin, int iflag)
{
   double tol(1.0e-4);
   double xmin, xmax;
   vector<double> LFparam;
   double lphi0, lL0, gamma1,gamma2;
   double phi0, L0;
   double alpha, beta;
   double scaled_integral, integral;
   LFparam=getLFParamQSO(z,iflag);
   lphi0=LFparam[0];
   lL0=LFparam[1];
   gamma1=LFparam[2];
   gamma2=LFparam[3];

   L0=pow(10.0,lL0);   // in Lsol
   phi0=pow(10.0,lphi0)/log(10.0);  // in Mpc^{-3}
   alpha= -(gamma1+1.0);
   beta= -(gamma2+1.0);
   cout<<"check: "<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<endl;

   //xmin=log(Lmin/L0);
   //xmax=log(1.0e3);   //stand in for infinity - should drop rapidly above L0 - CRAP!

   //cout<<"x: "<<xmin<<"\t"<<xmax<<endl;
//NEED TO REWRITE TO DO INTEGRAL IN dn/dLogL to get accurate results
//DOUBLE CHECK FACTORS OF L0 TOO

   //setDummyLFQSO_HOP(0.0,alpha,beta,this,1);
   //scaled_integral=qromb(dummyLFQSO_HOP,xmin,xmax,tol);

/////
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  //double xmin, xmax;
  ans = dvector(1,1);
  ans[1] = 0.0;
  setDummyLFQSO_HOP_ODE(L0,ans,ans,alpha,beta,this,1);
  xmin=Lmin/L0;
  cout<<"check lum: "<<Lmin<<"\t"<<L0<<endl;

  while (1) {
     step = xmin/10.0;
    xmax = xmin*10.0;
    oldans = ans[1];
    //cout<<"xx: "<<xmin<<"\t"<<xmax<<"\t"<<ans[1]<<"\t"<<oldans<<endl;
    odeint(ans,1,xmin,xmax,1.0e-4,step,0.0,&goodSteps,&badSteps,
	   dummyLFQSO_HOP_ODE,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    xmin = xmax;
  }
  scaled_integral = ans[1];
  free_dvector(ans,1,1);
/////

   
   integral=phi0*scaled_integral;
   return integral;
}

void setDummyLFQSO_HOP_ODE(double lL, double y[], double deriv[], double alphain, double betain, Galaxies *Gal1, int iflag)
{
   static double alpha, beta, L0;
   static Galaxies *Gal;
   double phi,L, phiband;
   double nu0(1.0);

   if(iflag==1){
      L0=lL;
      alpha=alphain;
      beta=betain;
      Gal=Gal1;
      return;
   }

   L=lL;
   phi=Gal->luminosityFunctionQSO_dNdL(L,1.0,1.0,alpha,beta);
   //convert bolometric phi to the B-band phi
   phiband=Gal->boloToBandPhi(nu0,L*L0*LSOLAR,phi);
   //cout<<L<<"\t"<<phi<<"\t"<<phiband<<endl;
   //cout<<L<<"\t"<<phi<<endl;
   phi=phiband;
   deriv[1]=phi;
}

void dummyLFQSO_HOP_ODE(double lL, double y[], double deriv[])
{
  setDummyLFQSO_HOP_ODE(lL,y,deriv,0.0,0.0,NULL,0);
}




double dummyLFQSO_HOP(double x)
{
   return setDummyLFQSO(x,0.0,0.0,NULL,0);
}

double setDummyLFQSO_HOP(double lL, double alphain, double betain, Galaxies *Gal1, int iflag)
{
   static double alpha, beta;
   static Galaxies *Gal;
   double phi,L;

   if(iflag==1){
      alpha=alphain;
      beta=betain;
      Gal=Gal1;
      return 0.0;
   }

   L=exp(lL);
   phi=L*Gal->luminosityFunctionQSO_dNdL(L,1.0,1.0,alpha,beta);
   return phi;
}



///////////////////////////////////////////////////////
// Hopkins et al QSO bolometric luminosity function
//////////////////////////////////////////////////////

//Calculate quasar luminosity function using a double power law form
// alternative form for d\Phi/d log L
// four free parameters for the fit:
//   phi0 : Mpc^{-3}
//   gamma1, gamma2 
//   L0   : ergs s^{-1}  or erg s^{-1} Hz^{-1}
// units: phi: no. galaxies per unit comoving volume per unit luminosity
double Galaxies::luminosityFunctionQSO_alt(double L, double phi0, double L0, double gamma1, double gamma2)
{
  double phi;
  double beta1, beta2;
  double phistar;

  beta1= -(gamma1+1.0);
  beta2= -(gamma2+1.0);
  phistar=phi0/log(10.0);
  phi=luminosityFunctionQSO(L,phistar,L0,beta1,beta2);
  return phi;
}

//Calculate quasar luminosity function using a double power law form
// for d\Phi/d L following Hopkins et al (2007) Equation (7)
//
// four free parameters for the fit:
//   phi0 : Mpc^{-3}
//   gamma1, gamma2 
//   L0   : ergs s^{-1}  or erg s^{-1} Hz^{-1}
// units: phi: no. galaxies per unit comoving volume per unit luminosity
double Galaxies::luminosityFunctionQSO_Hop(double L, double phi0, double L0, double gamma1, double gamma2)
{
  double phi;
  double alpha, beta;
  double phistarp;
  double LB, phiB;

  alpha= -(gamma1+1.0);  //convert to standard definition of first index
  beta= -(gamma2+1.0);  //convert to standard definition of second index
  phistarp=phi0/log(10.0); //renormalise phistar appropriately
  //phi=luminosityFunctionQSO(L,phistar,L0,beta1,beta2);
  phi=luminosityFunctionQSO_dNdL(L,phistarp,L0,alpha,beta);


  phiB=phi;

  //convert to bolometric luminosity
  //LB=boloToBandL(0.0,L);
  //phiB=boloToBandPhi(0.0,L,phi);

  return phiB;
}


/////////////////////////////////////////////////////////////////////////////
// Bolometric integrated luminosity to band luminosity conversion functions
// - based on Hopkins et al. (2007)
/////////////////////////////////////////////////////////////////////////////


// Apply correction factor for absorption column to bolometric
// luminosity function to get the band luminosity function
// Hopkins (2007) Eq (4)
//
// Units: L ergs s^-1
//        Lband erg s^{-1}
//        phi [Mpc^{-3}] is full luminosity function
double Galaxies::boloToBandPhi(double nu, double L, double phi)
{
   double f46, beta;
   double phiband;

   //B band conversion factors
   f46=0.260;
   beta=0.082;
   //B band is 4400 Angstrom

   phiband=phi*attenuationFactor(nu,L);

   return phiband;
}

// Apply correction factor for absorption column to bolometric
// luminosity function to get the band luminosity function
// Hopkins (2007) Eq (4)
//
// Units: L ergs s^-1
double Galaxies::attenuationFactor(double nu, double L)
{
   double f46, beta;
   double factor;

   // IR band conversion: 15 um
   f46=0.438;
   beta=0.068;

   // Soft X-ray: 0.5-2 keV
   f46=0.609;
   beta=0.063;

   //hard X-ray: 2-10 keV
   f46=1.243;
   beta=0.066;

   //B band conversion factors
   f46=0.260;
   beta=0.082;
   //B band is 4400 Angstrom

   factor=f46*pow(L/1.0e46,beta);

   return factor;
}


// Map bolometric L to band Luminosity Lband
// Hopkins (2007) Eq (2)
// Units: Lband ergs s^-1
//        L    ergs s^-1
//NOTE: currently only set up for B-band, but could be made
//      more general
double Galaxies::boloToBandL(double nu, double L)
{
   double c1, c2, k1, k2;
   double ratio, Lband;

   //HX band
   c1=10.83;
   k1= 0.28;
   c2=6.08;
   k2= -0.02;

   //SX band
   c1=17.87;
   k1= 0.28;
   c2=10.03;
   k2= -0.020;

   //IR band
   c1=7.40;
   k1= -0.37;
   c2=10.66;
   k2= -0.014;

   //B band
   c1=6.25;
   k1= -0.37;
   c2=9.00;
   k2= -0.012;

   ratio=c1*pow(L/1.0e10/LSOLAR,k1);
   ratio+=c2*pow(L/1.0e10/LSOLAR,k2);

   Lband=L/ratio;
   return Lband;
}


// Map band Luminosity Lband to bolometric luminosity L
// Inversion of Hopkins (2007) Eq (2)
// Also accounting for the attenuation Factor of Eq(4)
// Units: Lband ergs s^-1
//        L    ergs s^-1
double Galaxies::bandToBoloL(double nu, double Lband)
{
   double Lmin,Lmax,L;
   double tol(1.0e-5);

   //Bolomteric correction for B-band is <factor of 100
   //and bolometric L always larger than band L (by definition)
   Lmin=Lband;
   Lmax=Lband*1.0e5;

   if(boloToBandL(nu,Lmin)>Lband){
      cout<<"zriddr lower limit to large: Lmin="<<Lmin<<"\t"<<Lband<<endl<<endl;
      return 0.0;
   }else if(boloToBandL(nu,Lmax)<Lband){
      cout<<"zriddr upper limit too small: Lmax="<<Lmax<<"\t"<<boloToBandL(nu,Lmax)<<"\t"<<Lband<<endl;
      return 0.0;
   }

   //perform simple root finding for inversion
   L=zriddrConst(dummyBoloToBandL,Lband,Lmin,Lmax,tol);

   return L;
}

double setDummyBoloToBandL(double L,double nu1 ,Galaxies *Gal1,int iflag)
{
   double Lband;
   static double nu;
   static Galaxies *Gal;
   
   if(iflag==1){
      nu=nu1;
      Gal=Gal1;
      return 0.0;
   }

   Lband=Gal->boloToBandL(nu,L);
   Lband*=Gal->attenuationFactor(nu,L);

   return Lband;
}

//dummy wrapper for root finding on boloToBandL
double dummyBoloToBandL(double L)
{
   return setDummyBoloToBandL(L,0.0,NULL,0);
}


////////////////////////////////////////////////////////////////
// QSO LF parameters based on fits in the Hopkins + (2007)
////////////////////////////////////////////////////////////////

//Return best fit parameters for the QSO luminosity function
//use fit from Hopkins, Gordon, Hernquist (2007) - ApJ 654, 731
//
// This is a wrapper that allows for different methods to be used
//
// Return: [log10(phi0/Mpc^{-3}),Log10(L0/Lsol),gamma1,gamma2]
// where Lsol=3.9e33 erg s^{-1}
//
//Units: phi Mpc^-3
//
vector<double> Galaxies::getLFParamQSO(double z, int iflag)
{
   vector<double> param;

   if(iflag==0){
      param=getLFParamQSO_FULL(z);
   }else if(iflag==1){
      param=getLFParamQSO_PLE(z);
   }else if(iflag==2){
      param=getLFParamQSO_OBS(z);
   }else{
      cout<<"uninitialised LF param"<<endl;
   }
   return param;
}


//Return best fit parameters for the QSO luminosity function
//use fit from Hopkins, Gordon, Hernquist (2007) - ApJ 654, 731
//
// Use observational parameters from Table 2 of Hopkins (2007)
//
// Return: [log10(phi0/Mpc^{-3}),Log10(L0/Lsol),gamma1,gamma2]
// where Lsol=3.9e33 erg s^{-1}
//
//Units: phi Mpc^-3
//
vector<double> Galaxies::getLFParamQSO_OBS(double z)
{
   double lphi0, lL0, gamma1,gamma2;
   vector<double> param;

   //some specific values from Table 2 of Hopkins (2007)
   if(z<0.1){
      lphi0= -5.45;
      lL0=11.94;
      gamma1=0.868;
      gamma2=1.97;
   }else if(z<0.5){
      lphi0= -4.66;
      lL0=12.24;
      gamma1=0.6;
      gamma2=2.26;
   }else if(z<1.0){
      lphi0= -4.63;
      lL0=12.59;
      gamma1=0.412;
      gamma2=2.23;
   }

   //prepare vector of parameters
   param.push_back(lphi0);
   param.push_back(lL0);
   param.push_back(gamma1);
   param.push_back(gamma2);

   return param;
}

//Return best fit parameters for the QSO luminosity function
//use fit from Hopkins, Gordon, Hernquist (2007) - ApJ 654, 731
//
// Pure Luminosity Evolution from Eq. (9) & (10) + Table 3
   // doesn't give a great fit
//
// Return: [log10(phi0/Mpc^{-3}),Log10(L0/Lsol),gamma1,gamma2]
// where Lsol=3.9e33 erg s^{-1}
//
//Units: phi Mpc^-3
//
vector<double> Galaxies::getLFParamQSO_PLE(double z)
{
   double lphi0, lL0, gamma1,gamma2;
   vector<double> param;
   double kL1, kL2, kL3, zeta, lLstar0;
   double zref(2.0);

   lphi0= -4.733;
   lLstar0= 12.965;
   kL1=0.749;
   kL2= -8.03;
   kL3= -4.40;
   gamma1= 0.517;
   gamma2= 2.096;

   zeta=log10((1.0+z)/(1.0+zref));
   lL0=lLstar0+kL1*zeta+kL2*pow(zeta,2.0)+kL3*pow(zeta,3.0);


   //prepare vector of parameters
   param.push_back(lphi0);
   param.push_back(lL0);
   param.push_back(gamma1);
   param.push_back(gamma2);

   return param;
}


//Return best fit parameters for the QSO luminosity function
//use fit from Hopkins, Gordon, Hernquist (2007) - ApJ 654, 731
//
// Full Hopkins model from Eq (17)-(20)
//
// Return: [log10(phi0/Mpc^{-3}),Log10(L0/Lsol),gamma1,gamma2]
// where Lsol=3.9e33 erg s^{-1}
//
//Units: phi Mpc^-3
//
vector<double> Galaxies::getLFParamQSO_FULL(double z)
{
   double lphi0, lL0, gamma1,gamma2;
   vector<double> param;
   double kL1, kL2, kL3, zeta, lLstar0;
   double kg1, kg2_1, kg2_2;
   double gamma10, gamma20;
   double rhoBH0;
   double zref(2.0);

//parameters from Hopkins own code
//P0=-4.8250643; P1=13.035753;   P2=0.63150872; P3=-11.763560; P4=-14.249833; P5=0.41698725; P6=-0.62298947; P7=2.1744386;  P8=1.4599393;  P9=-0.79280099; P10=0.;P11=0.;P12=0.;P13=0.;P14=0.;

   lphi0= -4.825;
   lLstar0= 13.036;
   kL1=0.632;
   kL2= -11.76;
   kL3= -14.25;
   gamma10= 0.417;
   kg1= -0.623;
   gamma20= 2.174;
   kg2_1= 1.460;
   kg2_2= -0.793;
   rhoBH0=4.81;

   //PLE form
   zeta=log10((1.0+z)/(1.0+zref));

   lL0=lLstar0+kL1*zeta+kL2*pow(zeta,2.0)+kL3*pow(zeta,3.0);

   //additional evolution of faint and bright end
   gamma1=gamma10*pow(10.0,kg1*zeta);
   gamma2=gamma20*2.0/(pow(10.0,kg2_1*zeta)+pow(10.0,kg2_2*zeta));

   //WARNING: gamma2 shape does not agree with Hopkins paper!?!
   //too low at peak and comes down too rapidly at higher z

   //impose physical limit on gamma2 to prevent divergence
   if(gamma2<1.3) gamma2=1.3;

   //prepare vector of parameters
   param.push_back(lphi0);
   param.push_back(lL0);
   param.push_back(gamma1);
   param.push_back(gamma2);

   return param;
}

////////////////////////////////////////////////////////////////////
//Calculate the ionizing photon emissivity at the Lyman limit
// in usints of 10^24 erg s^-1 Mpc^-1
//
// Lmin = lower bound on luminosty
// nuObs is frequency at which luminosity function is measured
// phi0, L0, alpha describe the galaxy luminosity function of the sample
//
// This is epsilon_24 defined in Bolton & Haehnelt 2007
double Galaxies::emissivityQSO_Hop(double Lmin, double nuObs, double z, double phi0, double L0, double gamma1, double gamma2)
{
  double epsilon;
  double lamObs;
  double beta1,beta2;

  //epsilon=lumDensityLF_QSO_Hop(Lmin,phi0,L0,beta1,beta2);
  epsilon/=nuObs;

  //connect to 912 angstrom with a brokwn power law
  lamObs=SPEEDOFLIGHT_CGS/nuObs;
  epsilon*=pow(1050.0e-8/lamObs,-0.5);
  epsilon*=pow(912.0/1050.0,-1.5);

  epsilon/=1.0e24;
  return epsilon;
}


///////////////////////////////////////////////////////
// Stark, Loeb, Ellis (2008) galaxy luminosity function model
//////////////////////////////////////////////////////

//star formation rate
// mass : msol
// sfr : msol yr^{-1}
double Galaxies::sfrHalo(double mass, double z)
{
   double sfr;
   double fstar(0.13);   //best fit values from Section 4.1
   double duty(0.2);
   double tLT, tH;

   tH=2.0/(3.0*UNH*c->getH()*c->hubbleZ(z));
   tLT=duty*tH/YEAR;
   sfr=fstar*(c->getOmegab()/c->getOmega0())*mass/tLT;

   return sfr;
}

//conversion from galaxy sfr to luminosity 
//use relation of Madau+ 1998 (NEED TO CHECK DETAILS)
// assumes Salpeter IMF
//lum : ergs s^-1 Hz^-1  at 1500 Angstrom
// sfr : Msol yr^-1
double Galaxies::lumFromSFR(double sfr)
{
   double lum(8.0e27);

   lum*=sfr;

   return lum;
}

//inversion
double Galaxies::sfrFromLum(double lum)
{
   return lum/8.0e27;
}

//inversion to get mass from lum
double Galaxies::massFromLum(double lum, double z)
{
   double sfr, mass;

   sfr=sfrFromLum(lum);
   mass=sfr/sfrHalo(1.0,z);

   return mass;
}


//get luminosity function
double Galaxies::starkLF(double lum, double z)
{
   double mass, nL;

   mass=massFromLum(lum,z);
   nL=nm(mass,z,c);

   return nL;
}

///////////////////////////////////////////////////////
// Total amount of stars formed assuming LF and model for L(SFR)
//////////////////////////////////////////////////////

//Units: sfd  Msol Mpc^{-3}
double Galaxies::totalSFD(double zmin)
{
   double sfr, sfrold(0.0),sfd(0.0);
   double z,zmax(40.0);
   double dz(0.25);
   double dzdt;

   z=zmin;
   sfr=sfrDensity(z);

   while(z<zmax){
      z+=dz;
      sfrold=sfr;
      sfr=sfrDensity(z);
      dzdt=(1.0+z)*c->hubbleZ(z)*c->getH()*UNH*YEAR;
      sfd+=0.5*(sfr+sfrold)*dz/dzdt;
      //cout<<z<<"\t"<<sfr<<"\t"<<sfrold<<"\t"<<sfd<<endl;
   }

   return sfd;
}

//return the star formation rate density at a given redshfit
//Units: Msol Year^{-1} Mpc^{-1}
double Galaxies::sfrDensity(double z)
{
   double lfd;
   double Lmin, alpha, phi0, L0, M0, nuObs;
   vector<double> lfparam;
   double sfr;
   double mmin;
   int mflag(1);

      lfparam=getLFParam(z);
      M0=lfparam[0];
      phi0=lfparam[1];
      alpha=lfparam[2];
      nuObs=SPEEDOFLIGHT_CGS/(1900.0e-8);
      L0=lumFromAbsMag(M0,z)*nuObs;

      if(mflag==1){
      mmin=coolMass(c,z);
      }else if(mflag==2){
      mmin=coolMassH2(c,z);
      }else if(mflag==3){
      mmin=minIonMass(c,z);
      }  

      Lmin=lumFromSFR(sfrHalo(mmin,z))*nuObs;
      lfd=lumDensityLF(Lmin,alpha,phi0,L0);
      sfr=sfrFromLum(lfd)/nuObs;

      return sfr;
}

//Units MSol Mpc^{-3}
double Galaxies::reionSFD(double fesc)
{
   double sfd;

   sfd=1.7e6/fesc;

   return sfd;
}
