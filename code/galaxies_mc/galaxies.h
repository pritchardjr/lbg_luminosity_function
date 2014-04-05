// Galaxies.h
//
// Collate various useful galaxy calculations relating to reionization
//

#ifndef GALAXIES_H
#define GALAXIES_H


#include "dcosmology.h"
#include "astrophysics.h"
#include "vector.h"
#include "spline.h"
#include "spline2D.h"

const double LSOLAR(3.839e33);  //solar luminosity in ergs s^{-1}

class Galaxies {

 public:
  Galaxies();
  ~Galaxies();

  //initiation
  void initGalaxies(Cosmology *c1, Astrophysics *a1);
  void setCosmology(Cosmology *c1);
  //Member functions

  //LF
  double luminosityFunction(double L, double alpha, double phi0, double L0);
  double numberDensityLF(double Lmin, double alpha, double phi0, double L0);
  double lumDensityLF(double Lmin, double alpha, double phi0, double L0);
  double luminosityFunctionM(double M, double alpha, double phi0, double M0);
  double numberDensityLF_M(double Mmin, double alpha, double phi0, double M0);
  vector<double> getLFParam(double z);


  //magnitude conversions
  double magFromL(double L, double z);
  double lumFromMag(double m, double z);

  double absMagFromL(double L, double z);
  double lumFromAbsMag(double M, double z);
  double absMagFromM(double M, double z);

  //emissivity
  double emissivityGal(double Lmin, double nuObs, double z, double phi0, double L0, double alpha);
  double ndotGal(double e25, double alphaS, double fesc);
  double gammaGal(double e25, double alphaS, double fesc, double mfp, double z);
  double recombLevel(double z, double xi, double C);

  //quasars
  double emissivityQSO(double Lmin, double nuObs, double z, double phi0, double L0, double beta1, double beta2);
  double luminosityFunctionQSO(double L, double phi0, double L0, double beta1, double beta2);

  double luminosityFunctionQSO_dNdLogL(double L, double phistar, double L0, double gamma1, double gamma2);
  double luminosityFunctionQSO_dNdL(double L, double phistarp, double L0, double alpha, double beta);
  double luminosityFunctionQSO_dNdM(double M, double phistarpp, double M0, double alpha, double beta);
  double lumDensityLF_QSO(double Lmin, double phi0, double L0, double beta1, double beta2);
  double ndotQSO(double e24, double alphaS);
  double gammaQSO(double e24, double alphaS, double mfp, double z);

  //Hopkins LF fit to Luminosity-dependent density evolution
  double luminosityfunctionQSO_LDDE(double L, double z);

  //Hopkins Luminosity density
  double numberDensityQSO_HOP(double z, double Lmin, int iflag=0);


  //Hopkins LF
  double luminosityFunctionQSO_alt(double L, double phi0, double L0, double gamma1, double gamma2);
  double luminosityFunctionQSO_Hop(double L, double phi0, double L0, double gamma1, double gamma2);
  double boloToBandPhi(double nu, double Lbol, double phiBol);
  double attenuationFactor(double nu, double L);
  double boloToBandL(double nu, double L);
  double bandToBoloL(double nu, double Lband);

  //different fits to QSO LF parameters and wrapper
  vector<double> getLFParamQSO(double z, int iflag=0);
  vector<double> getLFParamQSO_OBS(double z);
  vector<double> getLFParamQSO_PLE(double z);
  vector<double> getLFParamQSO_FULL(double z);

  double emissivityQSO_Hop(double Lmin, double nuObs, double z, double phi0, double L0, double gamma1, double gamma2);

// Stark+ galaxy LF
  double sfrHalo(double mass, double z);
  double lumFromSFR(double sfr);
  double sfrFromLum(double lum);
  double massFromLum(double lum, double z);
  double starkLF(double lum, double z);

//total SFD
  double totalSFD(double zmin);
  double sfrDensity(double z);
  double reionSFD(double fesc);

 private:
  Cosmology *c;
  Astrophysics *a; 

};

//LF qso integration
double dummyLFQSO(double x);
double setDummyLFQSO(double L, double beta1in, double beta2in, Galaxies *Gal1, int iflag);

//LF HOP QSO integration
double dummyLFQSO_HOP(double x);
double setDummyLFQSO_HOP(double L, double alphain, double betain, Galaxies *Gal1, int iflag);

void dummyLFQSO_HOP_ODE(double lL, double y[], double deriv[]);
void setDummyLFQSO_HOP_ODE(double lL, double y[], double deriv[], double alphain, double betain, Galaxies *Gal1, int iflag);

//boloToBandL inversion
double setDummyBoloToBandL(double L,double nu1 ,Galaxies *Gal1,int iflag);
double dummyBoloToBandL(double L);


#endif

