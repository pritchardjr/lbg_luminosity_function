// astrophysics.h
// Contains prototypes and constants needed for calculating atomic
// levels for spectroscopy calculations

#ifndef ASTROPHYSICS_H
#define ASTROPHYSICS_H

#include "spline.h"

// Hated cgs units used for ease of connection to astronomy texts
const double BOHRA0(0.52918e-8);     // Bohr radius in cm 
//const double ELECTRONMASS(9.10938e-28);              //Electron mass in cgs
//const double PROTONMASS(1.6726e-24);              //Proton mass in cgs
const double PLANCKCONSTANT(6.62607e-27);     // h in cgs
const double RYDBERG_CGS(2.1798719e-11);    // Rydberg in cgs
const double RYDBERG(13.60569172);          // Rydberg in Electron volts
//const double SPEEDOFLIGHT_CGS(2.99792458e10);  //Speed of light in cgs
const double ELECTRONCHARGE_CGS(4.8033e-10);   //Electron charge in esu units
const double NULYMANLIMIT(3.28988e15);  //Lyman limit frequency in cgs
const double SIGMATHOMSON(6.65246e-25); //Thomson scattering x-section (cgs)
const double RADIATIONA(7.5657e-15); //BB energy density erg cm^-3 K^-4

const double YP(0.247);   //Primordial Helium mass fraction
const double FHE(YP/4.0/(1.0-YP));  //Primodial He/H ratio
const double SOLARMASS(1.989e33); //Mass of sun in grams
const double JANSKY(1.0e-23);  //erg s^-1 cm^-2 Hz^-1


//star formation parameters
//const double FSTAR(0.01);
//const double FESC(0.1);
//const double NION(30000.0);
const double CLUMPING(2.0);
const double ZSFR(45.0);  // must be less than 50.0
const double ZLOW(2.5);   // must be larger than 0.1
//const double ZREION(6.5); // redshift definitely reionized by
//const double ZSFR(24.89);


const int blflag(0);
const int lynmax(23); //Maximum level contributing to Lyn

const double XrayEmin(200.0);  //minimum x-ray energy
const double XrayEmax(1.0e4);  //maximum xray energy

const int MASSFCN(1);  //0: PS  1: ST 2: Jenkins

#include "dcosmology.h"

class Astrophysics {

 public:
  Astrophysics(Cosmology *c1, int popflag_in, int xrayflag_in, int lyaxray_in, double fxray_in);
  ~Astrophysics();

  //Member functions
  double getFSTAR(void);
  double getFESC(void);
  double getNION(void);
  double getNLYA(void);
  double getZeta(void);
  int getPopflag(void);
  int getLyaXray(void);
  double getFXRAY(void);
  double clumpingIntMEHR(double ldeltai);
  double volumeIntMEHR(double ldeltai);
  double getZReion();
  void setCosmology(Cosmology *c1);

  // Utility functions
  void initAstrophysics(double fstar, double fesc, double nion, double fx, double nlya, int popflag1, int xrayflag1, int lyaxray1, int sourceflag1, double starburst);
  void setFstar(double fstar);
  void setNION(double nion);
  void setFX(double fx);
  void setNLYA(double nlya);

  // Hydrogen Transition Probability function
  double transitionS(double n, double l, double j, double np, double lp,
		     double jp);
  double transitionA(double n, double l, double j, double np, double lp,
		     double jp);
  double transitionATotal(double n, double l, double j);
  double hydrogenFreq(double n, double np);
  double hydrogenWave(double n, double np);
  double hydrogenEnergy(double n, double np);
  double hydrogenMatrix(double n, double l, double np, double lp);
  double polyR(double n, double l, double np, double lp);
  double transitionP(double n, double l, double j, double np, double lp,
		     double jp);
  double transitionF(double n, double l, double j, double np, double lp,
		     double jp);

  //  Math functions for Hypergeometric Functions
  double intPochhammer(double a, double k);
  double term2F1(double a, double b, double c, double x);
  double factorial(double n);

  //Photoionisation cross sections
  double xsectionPhotoIonise(double nu, double Z);
  double gaunt1F(double x, double z);
  double mfpPhotoIonise(double energy, double z, double xfree);
  double xsectionPhotoIoniseGen(double nu, int z, int n);

  // Black body functions
  double uBlackBody(double nu, double T);
  double fluxBlackBody(double nu, double T);
  double nfluxBlackBody(double nu, double T);
  double ndensityBlackBody(double nu,double T);

  //Starburst functions
  double lumTotStarBurst(double SFR);
  double lumNuStarBurst(double nu,double SFR);
  double jNuStarBurst(double nu, double r, double z);
  double heatingNuStarBurst(double nu, double r, double z);
  double heatingStarBurst(double r, double z);

  //SNr functions
  double lumTotSNr(double SFR);
  double lumNuSNr(double nu,double SFR);

  //Mini-Quasar functions
  double lumTotMiniQuasar(double SFR);
  double lumNuMiniQuasar(double nu,double SFR);

  //Primary photoelectron fractions
  double fracPrimaryElectron(double x, int iflag);

  //lyman alpha flux
  double lyaFlux(double z);
  double sourceEmission(double nu, double z);
  double lymanZMax(double z, double n);
  double lyaRecycle(double n);
  double getZHII(double z);

  //xray emission 
  double xrayFlux(double z);
  double sourceEmissionXray(double E, double z);
  double xrayHeating(double z, double xe);
  double XrayTauSimp(double z, double zp, double E);
  double XrayTau(double z, double zp, double E);
  double lyafluxXray(double z, double xe);

  //radio flux
  double radioFlux(double nu, double z);
  double radioEmission(double nu, double z);

  //star formation
  double globalSFR(double z);
  double globalSFRsim(double z);
  double globalSFRcol(double z);

  //Compton heating rate
  double heatingCompton(double z, double Tgas, double xfree, Cosmology *c);
  double recombH(double Tgas);

  //spline fcoll
  double splineFColl(double zuse, Cosmology *c);
  double dfcdzFast(double z);

  //Spline free electron fraction history
  double splineXFree(double zuse, double z[], double xe[], int nmax, int iflag);
  double getXFree(double zuse);

  //Global temperature + ionization history
  void globalHistory(void);
  void getTIGM(double zin, double result[]);
  void getTIGMnum(double zin, double result[]);
  void getTIGMcont(double zin, double result[]);

  double findZReion(void);
  double findZReionDetail(void);

  double getTK(double z);
  double getXI(double z);
  double getXE(double z);

  //Calculate optical depth to SLS
  double getTauCMB();

  //Clumping factor
  double getClumpingMEHR(double Deltai, double z);
  double getVolumeMEHR(double Deltai, double z);
  void initPVParamMEHR(double z);
  double getMEHRindex(double z);
  double deltaIMEHR(double z, double xi);
  double clumpingMEHR(double z, double xi);

  //RECFAST related functions
  void callRECFAST(void);
  void useRECFAST(double z, double result[]);
  void getRECFAST(double z, double result[]);
  double splineRECFASTT(double zuse);
  double splineRECFASTX(double zuse);

  //halo formation functions
  double jeansMassFull(double z, double tk);
  double filterMassFull(double z);

  // HQB Radio loud AGN model
  double blackHoleMass(double mh, double z);
  double luminosityEddington(double mb);
  double luminosityRadioAGN(double nu, double mh, double z, double R);
  double numberRadioAGN(double R);
  double dutycycleRadioAGN(double z);
  double fluxRadioAGN(double nu, double mh, double z, double R);
  double getRfromFluxLimit(double nu, double mh, double z, double Flim);
  double fracBrighterR(double nu, double mh, double z, double Flim);
  double numberDensRadioAGNBrighterF(double nu, double z, double Flim);
  double numberRadioAGNBrighterF(double nu, double z, double Flim);
  double emissivityRadioAGN(double nu, double z);
  double numberCountsRadioAGN(double Flim, double zMin, double zMax, double nu);
  double differentialNumberCountsAGN(double S, double nu, double zobs);
  double dndsdvRadioAGN(double S, double nu, double z);
  double massCutRadioAGN(double mh, double z);

  //Number count code
  double fluxMoments(double sMin, double sMax, double zMin, double zMax, double nu, int order);
  double fluxMomentsDmDz(double z, double nu, double sMin1, double sMax, int order);
  double fluxMass(double mass,double z, double nu);
  double massFlux(double flux, double z, double nu);
  double haloFluxDensity(double z, double zSB, double M, double nu);
  double haloStarFormationRate(double mh);
  double lumRadioGalaxy(double nu, double sfr);
  double lumRadioFreeFree(double nu, double sfr);
  double differentialNumberCountsRG(double S, double nu, double zobs);


 protected:
  double omb;
  double h;
  double ombhh;

  double FSTAR;
  double FESC;
  double NION;
  double FXRAY;
  double NLYA;
  int popflag;

  int xrayflag;  //1: SB  2: SNR  3: mini-quasar
  int lyaxray; //use lya from xrays: 1: yes 0: no 2: use only lya from x-ray
  int sourceflag;

  double STARBURSTAGE;  //duration of starburst in years
  double ZREION;

  Cosmology *c;

  int globalHistoryFlag;
  Spline globalTK;
  Spline globalXI;
  Spline globalXE;
  Spline MEHRclumping;
  Spline MEHRvolume;

  //  int lyaHistoryFlag;
  //Spline globalLya;
  //Spline globalDJalphaXDZ;
};

//integrands for calculating heating rate
double getLambdaNuSB(double lE);
double setLambdaNuSB(double lE, double r1, double z1,Astrophysics *a1,int iflag);
void setDerivsTGas(double z, double tgas[],double dtgas[],Cosmology *c1,Astrophysics *a1,int iflag);
void derivsTGas(double z, double tgas[], double dtgas[]);

void setDerivsXFree(double z, double xe[],double dxe[],Cosmology *c1,Astrophysics *a1,int iflag);
void derivsXFree(double z, double xe[], double dxe[]);
void setDerivsHistory(double z, double y[],double dy[],Cosmology *c1,Astrophysics *a1,int iflag);
void derivsHistory(double z, double y[], double dy[]);
double setTauInt(double z, double zri, Cosmology *c1, Astrophysics *a1,int iflag);
double getTauInt(double z);

//clumping factor functions
double getPVNormMEHR(double C0);
double setPVNormMEHR(double beta1, double C0, double delta01,int iflag);
double PVMEHR(double Delta, double A, double beta, double C0, double delta0);
double setPVMEHR(double Delta, double A1, double beta1, double C01, double delta01, int iflag);
double getPVMEHR(double Delta);
double getDeltaPVMEHR(double Delta);
double getDelta2PVMEHR(double Delta);
double getlPVMEHR(double lDelta);
double getlDeltaPVMEHR(double lDelta);
double getlDelta2PVMEHR(double lDelta);
double getlDelta3PVMEHR(double lDelta);

double setDeltaIMEHR(double deltai, double xi1, Astrophysics *a1, int iflag);
double dummyDeltaIMEHR(double deltai);

//Lyman alpha flux
double setDJalphaDzDn(double n1, double z1, double zp, Cosmology *c1, Astrophysics *a1, int iflag);
double getDJalphaDzDn(double zp);
double setDJalphaDz(double z1,double zp,Cosmology *c1, Astrophysics *a1, int iflag);
double getDJalphaDz(double zp);

//Xray flux
double setDJxrayDzDE(double E, double z1, double zp, Cosmology *c1, Astrophysics *a1,int iflag);
double getDJxrayDzDE(double E, double zp);
double getDJxrayDz(double zp);
double getDJxrayDlE(double lE);

//xray heating
double getDLambdaXrayDz(double zp);
double getDLambdaXrayDlE(double lE);
double getDLambdaXrayDzDE(double E,double zp, Astrophysics *a1, int iflag);
//xray optical depth
double setXrayTau(double z1, double zp, double E1, Cosmology *c1, Astrophysics *a1,int iflag);
double getXrayTau(double zp);

//lya flux from xray excitation of HI
double getDJalphaXrayDz(double zp);
double getDJalphaXrayDlE(double lE);
double getDJalphaXrayDzDE(double E, double zp,Astrophysics *a1, Cosmology *c1, int iflag);

//radio flux
double setDJradioDz(double nu1, double z, double zp, Cosmology *c1, Astrophysics *a1, int iflag);
double getDJradioDz(double zp);

//radio AGN numbers
void setRadioNumberIntegrand(double m, double y[], double deriv[],double nu1, 
			     double Flim1, double z1, Cosmology *c1, Astrophysics *a1,int flag);
void radioNumberIntegrand(double m, double y[], double deriv[]);

double dummyNumberRadioAGN(double R);
double setDummyNumberRadioAGN(double R, Astrophysics *a1, int iflag);
double setNumberCountsAGNInt(double z, double nu1, double Flim1, Astrophysics *a1, int iflag);
double numberCountsAGNInt(double z);

//radio AGN emissivity
void setRadioEmissIntegrand(double m, double y[], double deriv[],double nu1, 
			     double z1, Cosmology *c1, Astrophysics *a1,int flag);
void radioEmissIntegrand(double m, double y[], double deriv[]);

//flux moment
void setFluxMomentIntegrand(double m, double y[], double deriv[],double nu1, 
			    double z1, int order,Cosmology *c1, Astrophysics *a1,int flag);
void fluxMomentIntegrand(double m, double y[], double deriv[]);
double setHaloFluxDensityInt(double z, double mh1, double nu1,Astrophysics *a1,Cosmology *c1, int iflag);
double haloFluxDensityInt(double z);

double setFluxMomentsInt(double z, double nu1, double sMin1, double sMax1, int order1, Astrophysics *a1, Cosmology *c1,int iflag);
double fluxMomentsInt(double z);
double setDiffNumberCountRGInt(double z, double nu1, double S1, Astrophysics *a1, Cosmology *c1, int iflag);
double diffNumberCountRGInt(double z);
double setDiffNumberCountAGNInt(double z, double nu1, double S1, Astrophysics *a1, Cosmology *c1, int iflag);
double diffNumberCountAGNInt(double z);
void setRadioAGNdNdSIntegrand(double m, double y[], double deriv[],double S1, double nu1, double z1, Cosmology *c1, Astrophysics *a1,int flag);
void radioAGNdNdSIntegrand(double m, double y[], double deriv[]);

//filter mass dummy function
double setFilterMassKernel(double z, Astrophysics *a1, int iflag);
double filterMassKernel(double z);

#endif
