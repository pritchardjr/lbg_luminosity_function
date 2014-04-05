// Miscellaneous programs for Lyman alpha forest
//
//


#include <math.h>
#include <iostream>
#include <fstream>
#include "astrophysics.h"
#include "dcosmology.h"
#include "dnumrecipes.h"
#include "spline.h"
#include "lymanforest.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////
// Constructor/Destructor 
///////////////////////////////////////////////////////////////////////////
Lyman::Lyman()
{
  lls_flag=1;

   temp_flag=0;

   double tau;
   tau=PI*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS/ELECTRONMASS/NU_LYA;
   tau*=FOSCILLATOR_LYA;
   //cout<<"sigmaLyman="<<tau<<endl;

   lambdaNLLS_flag=0;
   norm_NLLS=1.0;
   ci_flag=0;
   column_index=1.0;
}

Lyman::~Lyman()
{

}

///////////////////////////////////////////////////////////////////////////
// Initiation
///////////////////////////////////////////////////////////////////////////
void Lyman::initLyman(Cosmology *c1, Astrophysics *a1)
{

	a=a1;
	c=c1;

}

///////////////////////////////////////////////////////////////////////////
// Member functions
///////////////////////////////////////////////////////////////////////////


void Lyman::setLLSFlag(int lls_flag1)
{
  lls_flag=lls_flag1;
}

void Lyman::setTempDens(double T0, double beta)
{
   T0_in=T0;
   beta_in=beta;
   temp_flag=1;
}

// Problem with this is that it requires a value of gamma to set dNdz for LLS
//
void Lyman::setNormNLLS(double norm)
{
  norm_NLLS=norm;
}

void Lyman::setLambdaLLS(int LLS)
{
  lambdaNLLS_flag=LLS;
}

void Lyman::setColumnIndex(double index)
{
  ci_flag=1;
  column_index=index;

  if(index<-3.0){
    ci_flag=0;
    column_index=1.0;
  }

}

///////////////////////////////////////////////////////////////////////////
// Mean free path calculation given Gamma
///////////////////////////////////////////////////////////////////////////

// Calculate mfp given gamma_{-12}
//
// Units:  mfp 		Mpc
//         gamma	10^{-12} s^{-1}
double Lyman::mfpFromGamma(double z, double gamma)
{
	double mfp, deltai;

	if(z>6.5) return 0.0;

	deltai=deltaIonFromGamma(z,gamma);
	mfp=mfpFromDelta(z,deltai);

	return mfp;
}

// Calculate delta_i threshold for MHR00 model given gamma_12
//
//  Assumes clump is at 10^4 K
double Lyman::deltaIonFromGamma(double z, double gamma)
{
	double deltai;

	deltai=49.5*pow((1.0+z)/7.0,-3.0)*pow(gamma,2.0/3.0);

	return deltai;
}

double Lyman::mfpFromDelta(double z, double deltai)
{
  double mfp, fv;
  double fvmax(1.0-1.0e-7);
  fvmax=1.0-1.0e-7;
  static double zsave(1.0e30), deltai_save(1.0e30);
  static double lambda0;

  if(lambdaNLLS_flag==0){
    lambda0=6.0e6/(c->hubbleZ(z)*c->getH()*H0_CMSMPC)*(1.0+z);
    lambda0/=mfpCorrection(z);

  fv=a->getVolumeMEHR(deltai,z);
  if(fv>fvmax){  //calculation of fv gets shoddy when fv very close to 1

	mfp=lambda0/pow(1.0-fvmax,2.0/3.0);
	cout<<"warning: mfp>horizon: "<<mfp<<"\t"<<fv<<endl;
	return mfp;
	return 2000.0;
  }
  mfp=lambda0/pow(1.0-fv,2.0/3.0);
  }

  //cout<<pow(deltai/49.5*pow((1.0+z)/7.0,3.0),3.0/2.0)<<"\t";
  //cout<<mfp<<"\t";

  if(lambdaNLLS_flag==1){
    mfp=getLambda0(z,pow(deltai/49.5*pow((1.0+z)/7.0,3.0),3.0/2.0));
    mfp/=mfpCorrection(z);	
    mfp*=0.62;
  }

  //cout<<mfp<<endl;

  return mfp;
}

double Lyman::getLambda0(double z, double gamma)
{
  double lambda0;
  double Tgas(2.4e4); //recombining gas therefore at high T - matches FO05 
                      // little ambiuous what this should be
  //Tgas=1.1e4; //value that matches MHR00 mfp

    //connect to lyman alpha forest model
    lambda0=lambdaLLS(Tgas, gamma,z);

  return lambda0;
}

//taken from Furlanetto & Oh (2005) (A12)
double Lyman::mfpCorrection(double z)
{

  if(lls_flag==0) return 1.0;

  //need to be careful these criterea match that for
  // MEHR in astrophysics.cc

  if(z>5.99) return 2.67894;
  if(z>3.99) return 2.57162;
  if(z>2.99) return 2.0445;
  return 1.7274;	

}


///////////////////////////////////////////////////////////////////////////
// N_ion calculation given Gamma
///////////////////////////////////////////////////////////////////////////

//calculate emission rate of ionizing photons per unit comoving volume
// given gamma_{-12} 
// alphaS= source spectrum ; alphaB = reprocessed background spectrum
//
// Units:  Nion  s^{-1} Mpc^{-3}

double Lyman::nionFromGamma(double z, double gamma, double alphaS, double alphaB)
{
	double mfp, Nion;
	double colindex(1.0); //power law index of high N_{HI} absorbers

	//default spectrum choice
	if(alphaS<0.0 || alphaB<0.0){
		alphaS=3.0;
		alphaB=3.0;
	}

        a->initPVParamMEHR(z);
	mfp=mfpFromGamma(z,gamma);
	colindex=a->getMEHRindex(z)*2.0/3.0;
	if(ci_flag==1) colindex=column_index;
	Nion=pow(10.0,51.2)*gamma*(3.0/alphaS)*((alphaB+3.0*(2.0-colindex))/6.0);
	Nion*=(40.0/mfp)*pow((1.0+z)/7.0,-2.0);

	return Nion;
}

///////////////////////////////////////////////////////////////////////////
//  Inversion to get Gamma given Ndot
///////////////////////////////////////////////////////////////////////////
double Lyman::getGammaFromNdot(double z, double Ndot)
{
	double gamma;
	double gamma1(0.00001);
	double gamma2(50.0);
	double tol(1.0e-4);
	double nion1, nion2;

  if(z>6.5) return 0.0;

  a->initPVParamMEHR(z);
  setDummyGammaFromNdot(0.0,z,Ndot,this,1);

  //run into problems if gamma2 not large enough to bound nion
  nion1=dummyGammaFromNdot(gamma1);
  nion2=dummyGammaFromNdot(gamma2);
  if(fabs(nion1/fabs(nion1)-nion2/fabs(nion2))<1.0e-4){
    //cout<<"here "<<nion1<<"\t"<<nion2<<endl;
    while(fabs(nion1/fabs(nion1)-nion2/fabs(nion2))<1.0e-4){
      gamma2+=5.0;
      nion2=dummyGammaFromNdot(gamma2);
      //cout<<nion2<<endl;
    }
  }

  gamma=zriddrSimp(dummyGammaFromNdot,gamma1,gamma2,tol);

  return gamma;
}

double setDummyGammaFromNdot(double gamma, double z1, double Ndot1, Lyman *lyf1, int iflag)
{
	static Lyman *lyf;
	static double z, Ndot;
	double nion, mfp, deltai;

	if(iflag==1){
		lyf=lyf1;
		z=z1;
		Ndot=Ndot1;
		return 0.0;
	}

	nion= lyf->nionFromGamma(z,gamma,3.0,3.0);
	//mfp=lyf->mfpFromGamma(z,gamma);
	//deltai=lyf->deltaIonFromGamma(z,gamma);
	//cout<<gamma<<"\t"<<nion<<"\t"<<Ndot<<"\t"<<mfp<<"\t"<<deltai<<endl;

	return nion-Ndot;
}

double dummyGammaFromNdot(double gamma)
{
	return setDummyGammaFromNdot(gamma,0.0,0.0,NULL,0);
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Distribution of Lyman limit system
///////////////////////////////////////////////////////////////////////////
//  Calculate Column density of an object of overdensity \Delta=1+\delta
//  - associate overdensity to size of local Jeans Length
//  - assume ionization fraction set by photo-ionization equilibrium
//  - taken from Furlanetto and Oh (2005)
//
//Units:        NHI   cm^2

double Lyman::columnDensityFromDelta(double Delta, double Tgas, double gamma, double z)
{
   double NHI(3.3e17);

   NHI*=pow(Delta/100.0,1.5);
   NHI*=pow(Tgas/1.0e4,-0.26);     //Tgas in Kelvin
   NHI/=gamma;                   //gamma in 10^{-12} s^{-1}
   NHI*=pow((1.0+z)/7.0,9.0/2.0);

   return NHI;
}

double Lyman::deltaFromColumnDensity(double NHI, double Tgas, double gamma, double z)
{
   double Delta;

   Delta=columnDensityFromDelta(100.0,Tgas,gamma,z);
   Delta=100.0*pow(NHI/Delta,2.0/3.0);

   return Delta;
}

//ionization fraction in MHR00 model
double Lyman::xiFromDelta(double Delta, double Tgas, double gamma, double z)
{
   double xi;
   double chie(1.0+FHE);
   double gamma0(1.0e-12);

   xi=chie*c->nh(z)*recombRate(Tgas)/(gamma*gamma0)*Delta;
   //cout<<"xi="<<xi<<"\t"<<gamma<<endl;
   return xi;
}

//distribution of systems with column density N_{HI}
//This is defined so that
// dz/dl \int NHI f(NHI,z) dNHI = x_HI n_H
double Lyman::fColumnDensity(double NHI, double Tgas, double gamma, double z)
{
   double Delta, fNHI, xHI,dNHIdDelta;
   double fudge;

   Delta=deltaFromColumnDensity(NHI,Tgas,gamma,z);
   dNHIdDelta=1.5*NHI/Delta;
   xHI=xiFromDelta(Delta,Tgas,gamma,z);  //really xHI=1 for LLS
   //if(NHI>1.0/SIGMAION) xHI=1.0;
   //if(NHI>1.0/SIGMAION) dNHIdDelta=0.5*NHI/Delta;

   fNHI=c->getOmbhh()*CRITDENSITY*(1.0-YP)/PROTONMASS*SPEEDOFLIGHT_CGS;
   fNHI/=UNH*c->getH();
   fNHI/=dNHIdDelta;
   fNHI*=Delta*getPVMEHR(Delta);
   fNHI*=xHI/NHI;

   return fNHI*norm_NLLS;
}

//Number density of objects of column density NHI per unit NHI and z
// i.e. d^2N/dNHIdz
double Lyman::dNColumnDensity(double NHI, double Tgas, double gamma, double z)
{
   double dN;

   dN=fColumnDensity(NHI,Tgas,gamma,z);
   dN*=pow(1.0+z,2.0)/c->hubbleZ(z);
   return dN;
}

//Number density of LLS
double  Lyman::dNLLSdz(double Tgas, double gamma, double z)
{
   double dN;
   double lDMin;
   double lDMax;
   double tol(1.0e-5);
   double NHILLS(1.0/SIGMAION);


   lDMin=log(deltaFromColumnDensity(NHILLS,Tgas,gamma,z));
   lDMax=log(exp(3.5*log(10.0)));

   //cout<<exp(lDMin*log(10.0))<<endl;

   setDummyDNLLS(0.0,Tgas,gamma,z,this,1);

   dN=qromb(dummyDNLLS,lDMin,lDMax,tol);

   return dN;
}

double setDummyDNLLS(double lDelta, double Tgas1, double gamma1, double z1, Lyman *lyf1, int iflag)
{
   static double Tgas, gamma,z;
   static Lyman *lyf;
   double NHI, dN, Delta, dNHIdDelta;

   if(iflag==1){
      Tgas=Tgas1;
      z=z1;
      gamma=gamma1;
      lyf=lyf1;
      return 0.0;
   }

   Delta=exp(lDelta);
   //Tgas=lyf->tempDensityRelation(Delta,z);
   NHI=lyf->columnDensityFromDelta(Delta,Tgas,gamma,z);
   dN=lyf->dNColumnDensity(NHI,Tgas,gamma,z);
   dNHIdDelta=1.5*NHI/Delta;

   //cout<<Delta<<"\t"<<NHI<<"\t"<<dN<<"\t"<<dN*dNHIdDelta<<"\t"<<dN*dNHIdDelta*Delta<<endl;

   return dN*dNHIdDelta*Delta;
}

double dummyDNLLS(double lNHI)
{
   return setDummyDNLLS(lNHI,0.0,0.0,0.0,NULL,0);
}

double Lyman::lambdaLLS(double Tgas, double gamma, double z)
{
  double lamLLS, dNdz;

   lamLLS=SPEEDOFLIGHT_CGS;
   lamLLS/=c->hubbleZ(z)*UNH*c->getH();
   dNdz=dNLLSdz(Tgas,gamma,z);
   lamLLS/=dNdz;
   //lamLLS/=dNLLSdz(Tgas,gamma,z);
   //cout<<dNLLSdz(Tgas,gamma,z)<<endl;
   //cout<<"gamma= "<<gamma<<"\t"<<"dNdz= "<<dNdz<<"\t"<<"lamdaLLS= "<<lamLLS/MPC<<endl;
   return lamLLS/MPC;
}

///////////////////////////////////////////////////////////////////////////
//  Conversion of tau_eff to gamma_ion  from FG08
///////////////////////////////////////////////////////////////////////////

//optical depth for lyman alpha photons
//assumes fully ionized medium
//Delta=1+delta
double Lyman::tauLya(double z, double T, double gamma, double Delta)
{
   double tau;
   double gamma0(1.0e-12);
   double zHEII(3.0);
   

   tau=PI*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS/ELECTRONMASS/NU_LYA;
   tau*=FOSCILLATOR_LYA;

   tau/=c->hubbleZ(z)*UNH*c->getH();
   tau*=recombRate(T)*c->nh(z)*c->nh(z)/(gamma*gamma0);
   tau*=Delta*Delta;
   if(z<zHEII){
      tau*=1.0+2.0*FHE;
   }else{
      tau*=1.0+FHE;
   }

   return tau;
}

//hydrogen recombination rate
double Lyman::recombRate(double T)
{
   double T0(1.0e4);            // Kelvin
   //double alphaA(1.0);
   double alphaA(4.2e-13);   // cm^3 s^-1
   double alphaB(1.6e-13);

   return alphaA*pow(T/T0,-0.7);
   return alphaB*pow(T/T0,-0.7);
   
}

//Hui & Gnedin temperature density relation
//Delta=1+delta
//T0=2.2/pm0.2 *10^4 K from Zaldarriaga et al 2001
//beta=0.62 for early reionization
//note gamma=beta+1 for adiabatic index gamma
double Lyman::tempDensityRelation(double Delta, double z)
{
   double T;
   double beta(0.62);
   double T0(2.2e4);  //Kelvin
   double zm, zp, Tm, Tp;

   if(z<=3.0){
      Tp=2.3e4;
      Tm=2.1e4;
      zp=3.0;
      zm=2.4;

      T0=Tm+(z-zm)*(Tp-Tm)/(zp-zm);
   }else if(z<=3.9){
      Tp=2.2e4;
      Tm=2.3e4;
      zp=3.9;
      zm=3.0;

      T0=Tm+(z-zm)*(Tp-Tm)/(zp-zm);      
   }else{
      T0*=pow((1.0+z)/4.9,2.0);
   }

   if(BOLTON_FLAG==1){
      T0=1.0e4;
      beta=0.3;
   }

   if(temp_flag==1){
      T0=T0_in;
      beta=beta_in;
   }

   T=T0*pow(Delta,beta);

   return T;
}

//Mean transmittance of flux
double Lyman::meanTransmittance(double z, double gamma)
{
   double meanF;
   double lDeltamin(-14.0);
   double lDeltamax(12.0);
   double tol(1.0e-4);

   setMeanFKernel(z,gamma,this,1);
   //a->initPVParamMEHR(z);
   meanF=qromb(getMeanFKernel,lDeltamin,lDeltamax,tol);

   return meanF;
}

double setMeanFKernel(double lDelta, double gamma1, Lyman *lyf1, int iflag)
{
   static Lyman *lyf;
   static double z;
   static double gamma;
   double tau,T,Delta;

   if(iflag==1){
      z=lDelta;
      lyf=lyf1;
      gamma=gamma1;
      return 0.0;
   }

   Delta=exp(lDelta);
   T=lyf->tempDensityRelation(Delta,z);
   tau=lyf->tauLya(z,T,gamma,Delta);
   //cout<<Delta<<"\t"<<T<<"\t"<<tau<<endl;
   return exp(-tau)*getPVMEHR(Delta)*Delta;
}

double getMeanFKernel(double lDelta)
{
   return setMeanFKernel(lDelta,0.0,NULL,0);
}

double Lyman::tauEff(double z, double gamma)
{
   return -log(meanTransmittance(z,gamma));
}

//convert observed values of effective tau
// into value for gamma_{-12}
double Lyman::gammaFromTau(double z, double tau)
{
   double gamma;
   double gamma1(0.01);
   double gamma2(5.0);
   double tol(1.0e-4);

   a->initPVParamMEHR(z);
   setGammaFromTau(tau,z,this,1);

   //cout<<gamma1<<"\t"<<gamma2<<"\t"<<getGammaFromTau(gamma1)<<"\t"<<getGammaFromTau(gamma2)<<endl;
   if(getGammaFromTau(gamma1)<0.0){
     while(getGammaFromTau(gamma1)<0.0){
     gamma1/=1.05;
     //cout<<"reducing gamma1: "<<gamma1<<"\t"<<getGammaFromTau(gamma1)<<endl;
     }
   }
   gamma=zriddrSimp(getGammaFromTau,gamma1,gamma2,tol);
   return gamma;
}

double setGammaFromTau(double gamma, double z1, Lyman *lyf1, int iflag)
{
   static double z;
   static Lyman *lyf;
   static double tauIn;
   double tau;

   if(iflag==1){
      tauIn=gamma;
      z=z1;
      lyf=lyf1;
      return 0.0;
   }

   tau=lyf->tauEff(z,gamma);
   //cout<<gamma<<"\t"<<tau<<"\t"<<tauIn<<endl;
   return tau-tauIn;
}

double getGammaFromTau(double gamma)
{
   return setGammaFromTau(gamma,0.0,NULL,0);
}
