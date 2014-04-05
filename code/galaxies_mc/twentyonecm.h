// twentyonecm.h
// Contains prototypes and constants needed for calculating 21cm signal

#ifndef TWENTYONECM_H
#define TWENTYONECM_H

const double FALPHA(0.4162);         //lyman alpha oscillator strength
const double A21CM(2.85e-15);
const double nu21cm(1.420405e9);
const double lambda21cm(21.12);
const double T21cm(nu21cm*PLANCKCONSTANT/BOLTZK);

const double KCUT(500.0);
//const double KCUT(100.0);

#include "dcosmology.h"
#include "astrophysics.h"
#include <vector>
#include "reionization.h"

class TwentyOneCM {

 public:
  TwentyOneCM(Cosmology *c1, Astrophysics *a1);
  ~TwentyOneCM();

  //Member functions


  //temperatures
  double tBright(double z);
  double tSpin(double z);
  double tau21CM(double z, double tspin, double xi);
  double tKinetic(double z);
  double getXFree(double z);
  double tBrightGen(double z, double tk, double xi, double lya);
  double tSpinGen(double z, double tK, double Xi, double lya);
  double tBrightSat(double z, double xi);
  double tSky(double z);

  //coupling
  double getXColl(double z, double tk, double xi);
  double getXCollEH(double z, double tk, double xi);
  double getXCollHH(double z, double tk, double xi);
  double getXAlpha(double z, double tk, double xi, double lyaflux);
  double getSAlpha(double z, double tk, double xi);

  //fluctuation coefficients
  void getBeta(double z, double tk, double xi, double lyaflux, double beta[]);
  void getBetaTb(double z, double tk, double xi, double lyaflux, double betaT[]);

  //Adiabatic index
  double getGammaA(double z, double tk);

  //lya poisson fluctuation code
  double poissonPowerLya(double z);
  double weightPLya(double z,double zp, double tk, double xi);
  double poissonCorrelationLya(double l, double z);

  //collision kappas
  int Zygelman_init(void);
  double kappaHH(double T);
  double kappaEH(double tk);
  double getDLogKappaHH(double tkin);
  double getDLogKappaEH(double tkin);

  //21cm cutoffs
  double getCutoffThermal(double z, double tk);
  double getCutoffFilter(double z);

  //Lya fluctuation code
  double getWindowLya(double z, double k);
  double powerSpectrumLya(double z, double k, double store[]);

  //halo averaged bias
  double splineBias(double z);
  double meanBiasPS(double z);
  double meanBiasST(double z);
  double splineMFColl(double xuse);


  //xray fluctuations
  double powerSpectrumXray(double z, double k, double store[]);
  double getWindowXray(double z, double k);

  //xray poisson fluctuation code
  double poissonPowerXray(double z);
  double weightPXray(double z,double zp, double tk, double xi);
  double poissonCorrelationXray(double l, double z);

  //T fluctuation evolution
  double QXray(double z);
  double QCompton(double z);
  double QIon(double z);
  double QRec(double z);
  double getGTHistory(double zin,double k);
  double getGT(double zin,double k, double gT[]);
  double splineWkXray(double zuse, double k, int useflag);
  void getGTN(double zin[], double k, double **gT, int nz);

  //Brightness temperature power spectrum
  double powerSpectrumFull(double z, double k, double store[], int iflag);
  void powerSpectrumFullN(double zin[], double k, double **store,int nz, int iflag);
  
  //Kinetic temperature power spectrum
  double powerSpectrumTemp(double z, double k, double store[]);
  void powerSpectrumTempN(double zin[], double k, double **store, int nz);

  //Radio fluctuations
  double getWindowRadio(double z, double k);
  double getBetaTbRadio(double z, double tk, double xi, double lyaflux);
  double getBetaTbRadioBack(double z, double tk, double xi, double lyaflux);

  //post reionization
  double postReionBias21CM(double z, double R);
  double postReionPowerSpectrum(double z, double R, double powerDD);

  //calculating full power spectrum
  double fullPowerSpectrum(double z, double k, double powerFZH, double powerDD, double powerXD, double betaD, double betaX, double betaL, double betaT, double betaV, double gt, double wa, double tk, double xi, double rF, double rT, double powerPostTB);

  //transtion redshift for astrophysics to physics
  double transitionZMeanLya();
  double transitionZMeanT();
  double transitionZLya(double k);
  double transitionZT(double k);

// postreionization code
	double tBrightPR(double z, double delta, double R);
	double fillingMass(double deltai, double delta, double R);

// 21 cm pdf calculation from FZH04
        double getPsiPDF(double psi, double z, double zeta, double mpix);
        double XIfromOverdensity(double delta, double z, double zeta, double mpix);
        double psiFromDelta(double delta, double z, double zeta, double mpix);
        double derivPsiByDelta(double delta, double z, double zeta, double mpix);
        vector<double> deltaFromPsi(double psi, double z, double zeta, double mpix);
        double probDelta(double delta, double z, double mpix);
        double getBrightnessPDF(double psi, double z, double zeta, double mpix);
        double getProbPixExternalIon(double sigm, double delta, double z, double zeta, double mpix);
        double evaluateCondProb(double delta, double z, double zeta, double mpix);


 protected:
  Cosmology *c;
  Astrophysics *a;
};


double setWindowLyaK(double z1, double zuse, double k1,double n1,double tk1,Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag);
double intWindowLyaK(double zuse);

double setWindowXrayK(double z1, double zuse, double k1,double tk1,Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag);
double intWindowXrayK(double zuse);
void setDerivsGT(double z, double k1, double y[],double dy[],Cosmology *c1,TwentyOneCM *tocm1,int iflag);
void derivsGT(double z, double y[], double dy[]);

//poisson code -Lya
double setPoissonCorrelationKernelLya(double rA, double theta, double z1,double l1, double tk1, Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag);
double intPoissonLya(double rA, double theta);
double intPoissonRLya(double lrA);
double intPoissonThetaLya(double theta);

double fMColl(double z, double mMin, Cosmology *c,int massFcn);
double setFMCollInt(double mass, Cosmology *c1, double z1, int fcn, int flag);
double fMCollInt(double mass);

void FTTable(double xivec[], double rvec[], int nr1, double pvec[], 
	     double kvec[], int nk1);
double setFTQRInt(double r, double kp, int flag);
double FTQRInt(double r);

//poisson code - Xray
double setPoissonCorrelationKernelXray(double rA, double theta, double z1,double l1, double tk1, Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag);
double intPoissonXray(double rA, double theta);
double intPoissonRXray(double lrA);
double intPoissonThetaXray(double theta);

double setWindowRadioK(double z1, double zuse, double k1,Cosmology *c1, Astrophysics *a1, TwentyOneCM *tocm1, int iflag);
double intWindowRadioK(double zuse);

//transition redshift code
double setDummyTransitionZMeanLya(double z, Astrophysics *a1, TwentyOneCM *tocm1, int iflag);
double dummyTransitionZMeanLya(double z);
double setDummyTransitionZMeanT(double z, Astrophysics *a1, TwentyOneCM *tocm1, int iflag);
double dummyTransitionZMeanT(double z);
double setDummyTransitionZLya(double z, Spline *wLya1, Spline *wDens1, int iflag);
double dummyTransitionZLya(double z);
double setDummyTransitionZT(double z, Spline *wTemp1, Spline *wDens1, int iflag);
double dummyTransitionZT(double z);

//post reionization code

double setFillingMassInt(double z, double deltai1, double delta1, double zobs1, Cosmology *c1, Astrophysics *a1, int iflag);
double fillingMassInt(double z);

//pdf code
double setGetPsi(double delta, double z1, double zeta1, double mpix1, TwentyOneCM *tocm1, int iflag);
double dummyGetPsi(double delta);

double getCondProbKernel(double sigm);
double setCondProbKernel(double sigm, double delta, double z, double zeta, double mpix, TwentyOneCM *tocm1, int iflag);

#endif
