//lymanforest.h
//misc functions relating to the lyman lapha forest

#ifndef LYMAN_H
#define LYMAN_H

const double FOSCILLATOR_LYA(0.4162);
const double NU_LYA(2.46741e15);        //Lya frequency in second^-1
const double SIGMAION(6.3e-18); //HI Ionization cross-section at Lyman limit cm^2
const double SIGMALYMAN(7.03e-11); //Lyman alpha scattering cross-section on resonance in cgs
const int BOLTON_FLAG(1);    //use BH07 values in temperature density


class Lyman {

 public:
  Lyman();
  ~Lyman();

 //initiation
  void initLyman(Cosmology *c1, Astrophysics *a1);

  //member functions
  void setLLSFlag(int lls_flag1);
  void setTempDens(double T0, double beta);
  void setNormNLLS(double norm);
  void setLambdaLLS(int LLS);
  void setColumnIndex(double index);

 //mfp calculation
  double mfpFromGamma(double z, double gamma);
  double deltaIonFromGamma(double z, double gamma);
  double mfpFromDelta(double z, double deltai);
  double getLambda0(double z, double gamma);
  double mfpCorrection(double z);

 //Nion calculation
  double nionFromGamma(double z, double gamma, double alphaS, double alphaB);
  double getGammaFromNdot(double z, double Ndot);

  //Lyman limit system distribution
  double columnDensityFromDelta(double Delta, double Tgas, double gamma, double z);
  double deltaFromColumnDensity(double NHI, double Tgas, double gamma, double z);
  double xiFromDelta(double Delta, double Tgas, double gamma, double z);
  double fColumnDensity(double NHI, double Tgas, double gamma, double z);
  double dNColumnDensity(double NHI, double Tgas, double gamma, double z);
  double dNLLSdz(double Tgas, double gamma, double z);
  double lambdaLLS(double Tgas, double gamma, double z);
 
//conversion from tau_effective to gamma_ion
  double tauLya(double z, double T, double gamma, double Delta);
  double recombRate(double T);
  double tempDensityRelation(double Delta, double z);
  double meanTransmittance(double z, double gamma);
  double tauEff(double z, double gamma);
  double gammaFromTau(double z, double tau);

 protected:

  int lls_flag;

  Cosmology *c;
  Astrophysics *a;

  int temp_flag;  //fix temperature-density relation
  double T0_in;
  double beta_in;

  double lambdaNLLS_flag;  //use dNdz for LLS to set mfp
  double norm_NLLS;  //normalise dNLLSdz to data

  int ci_flag; //let column index be set externally
  double column_index;

  
};

//LLS distribution
double setDummyDNLLS(double lDelta, double Tgas1, double gamma1, double z1, Lyman *lyf1, int iflag);
double dummyDNLLS(double lDelta);

//invert Nion to get gamma
double setDummyGammaFromNdot(double gamma, double z1, double Ndot1, Lyman *lyf1, int iflag);
double dummyGammaFromNdot(double gamma);


//tau_effective
double setMeanFKernel(double Delta, double gamma1, Lyman *lyf1, int iflag);
double getMeanFKernel(double Delta);
double setGammaFromTau(double gamma, double z1, Lyman *lyf1, int iflag);
double getGammaFromTau(double gamma);
#endif
