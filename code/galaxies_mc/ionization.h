// Ionization.h
//
// Simplified reionization code for use in project with Avi and Stuart
//

#ifndef IONIZATION_H
#define IONIZATION_H


#include "dcosmology.h"
#include "astrophysics.h"
#include "vector.h"
#include "spline.h"
#include "lymanforest.h"
#include "spline2D.h"

const double ZION_START(40.0);
//const double ZION_END(4.0);
const double ZION_END(3.0);
const double clumping_flag(1);  //use MHR clumping model

class Ionization {

 public:
  Ionization(Cosmology *c1, Astrophysics *a1);
  ~Ionization();

  //Member functions
  void setFlags(int ndot, int cmb, int lls, int nt);
  double getZeta(double z);  //dummy function for Zeta
  double getZCut();
  double getClumpingParam();
  double setZetaParam(double z, vector<double> zparam, int iflag);
  double getNion(double z);
  void setPDF(string file);
  void setGammaFlag(int gamma_flag1);
  void setTocmFlag(int tocm1);
  void setNocutFlag(int nocut1);
  void setNocmbFlag(int nocmb1);

  //polynomial code
  double getGamma4();
  double checkNionMax(double Nion, double a, double b, double c, double zcut);

  //ionization
  double getXI(double zin);
  void globalHistory(void);
  double findZReion();

  //tau
  double getTauCMB();

  //likelihood
  double likelihoodZeta(vector<double> param);
  double gaussianProb(double x, double x0, double sig);
  double uniformProb(double x, double xmin, double xmax);
  double nongaussianProb(double x, Spline *spline, string file, int iflag);

  //clumping
  void setClumping();
  double getClumping(double z, double xi);

 private:
  Cosmology *c;
  Astrophysics *a; 
  double zetaBub;

  double zion_cut;
  double clumping_param;
  double ZREION;

  Spline globalXI;
  Spline2D clumpingSP;
  int globalHistoryFlag;

  int cmb_flag;
  int ndot_flag;
  int lls_flag;
  int nthread;
  string filetag;

  //splines for Ndot likelihoods
  Spline ndot4sp;
  Spline ndot5sp;
  Spline ndot6sp;
  string pdfdir;
  int gamma_flag;
  int tocm_flag;
  int nocut_flag;
  int nocmb_flag;

};

//polynomial 
double polynomial(double x, double a, double b, double c);
double checkPolynomial(double x, double a, double b, double c);
double dummyGamma(double gamma);
double setDummyGamma(double gamma, double zcut1, double a1, double b1, int iflag);

void setDerivsHistoryIon(double z, double y[],double dy[],Cosmology *c1,Astrophysics *a1, Ionization *Ion1, int iflag);
void derivsHistoryIon(double z, double y[], double dy[]);

//tau
double setTauIntIon(double z, double zri1,Cosmology *c1, Ionization *Ion1,int iflag);
double getTauIntIon(double z);

//vegas integrand
double vegasIntegrand(double param[], double wt);
double setVegasIntegrand(double paramin[], Ionization *Ion1, int param_flag1, string file, int iflag);
#endif
