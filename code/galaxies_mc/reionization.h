// reionization.h
// Contains prototypes and constants needed for calculating HII regions
// during reionization
//
//
// References: FZH04 - Furlanetto, Zaldarriaga, Hernquist astro-ph/0403697
//

#ifndef REIONIZATION_H
#define REIONIZATION_H


#include "dcosmology.h"
#include "astrophysics.h"

//const double ZETABUB(45.0);

class Reionization {

 public:
  Reionization(Cosmology *c1, Astrophysics *a1);
  ~Reionization();

  //Member functions
  double getZeta();
  void setZeta(double zetaIn);
  double getXICUT();
  
  //fudge function
  void normZeta(double z);

  //Bubble distribution code
  double nmBub(double mb, double z);
  double fillingQ(double z);
  double barrierBub(double mb, double z, double BB[]);
  double meanRBub(double z);
  double meanNBub(double z);
  double biasBub(double mb, double z);
  double biasBubFZH(double mb, double z);
  double meanBiasBub(double z);

  //Correlation function of HII regions
  double correlateTb(double r, double z);
  double correlateXX(double r, double z);
  double correlateXD(double r, double z);
  double correlateMeanBub(double r, double mb, double z);
  double overlapVolume(double r, double mb);
  double meanBubOverlap(double r, double z);
  double meanBubOverlapBias(double r, double z);
  double meanBubOverlapCluster(double r, double z);
  double correlateXX1b(double r, double z);
  double correlateXX2b(double r, double z);
  double correlateXDIntegral(double r, double z);
  double meanXiXD(double r, double z);
  double xiAverageInt(double mb, double r, double z);
  double xiAverageInt2(double mb, double r, double z);

  //McQuinn related code
  double correlateXXMQ(double r, double z);
  double correlateXDMQ(double r, double z);
  double p2MQ(double r, double z);

  // power spectrum code
  int getRVec(double rvec[]);
  double getCorrDD(double z, double rvec[], double xi[], int nr);
  double getCorrXD(double z, double rvec[], double xi[], int nr);
  double getCorrXX(double z, double rvec[], double xi[], int nr);
  double getCorrTB(double z, double rvec[], double xi[], int nr, double corrDDV[], double corrXDV[], double corrXXV[], double xH);

  //simple FZH power spectrum
  double getPowerFZH(double z, double k, double xi);

  //clumping factor code
  double clumpingBub(double z, double xi);
  double deltaIBub(double mb, double z);

 private:
  Cosmology *c;
  Astrophysics *a;
  double zetaBub;
  double XICUT;
  int ASTROPHYSICS_FLAG;  //used to normalise bubbles to astrophysics history
};

//inverse error function code
double inverseERF(double y);
double findIERF(double x);
double setIERF(double x, double y1, int iflag);

//

void setQIntegrand(double m, double y[], double deriv[], 
		   double z1, Reionization *rz1, int flag);
void qIntegrand(double m, double y[], double deriv[]);

double setRBub(double lmb, double z1, Reionization *rz1, Cosmology *c1, int iflag);
double intRBub(double lmb);
double intNBub(double lmb);

void setBiasBubIntegrand(double m, double y[], double deriv[], 
			 double z1, Reionization *rz1, int flag);
void biasBubIntegrand(double m, double y[], double deriv[]);

void setBubOverlapIntegrand(double m, double y[], double deriv[], 
			    double z1, double r1, Reionization *rz1, int flag);
void bubOverlapIntegrand(double m, double y[], double deriv[]);
void setBubOverlapClusterIntegrand(double m, double y[], double deriv[], 
			    double z1, double r1, Cosmology *c1, Reionization *rz1, int flag);
void bubOverlapClusterIntegrand(double m, double y[], double deriv[]);

void setCorrXDIntegrand(double m, double y[], double deriv[], 
			double z1, double r1, Cosmology *c1, Reionization *rz1, int flag);
void corrXDIntegrand(double m, double y[], double deriv[]);
void setXiXDIntegrand(double m, double y[], double deriv[], 
		      double z1, double r1, Cosmology *c1, Reionization *rz1, int flag);
void xiXDIntegrand(double m, double y[], double deriv[]);

double setXiAverageIntegrand(double rs, double mu, double r1, double z1, Cosmology *c1,int iflag);
double intXiAverageMu(double mu);
double intXiAverageR(double rs);
double intXiAverageMu2(double mu);
double intXiAverageR2(double rs);

//McQuinn related code
void setBubOverlapBiasIntegrand(double m, double y[], double deriv[], 
				double z1, double r1, Cosmology *c1,Reionization *rz1, int flag);
void bubOverlapBiasIntegrand(double m, double y[], double deriv[]);



//My 3D FT code
void FTTableJP(double xivec[], double rvec[], int nr1, double pvec[], 
	       double kvec[], int nk1);
double setFTQRIntJP(double r, double kp, double rmin1, double rmax1,int flag);
double FTQRIntJP(double r, void * params);

//clumping factor code
void setClumpingBubIntegrand(double m, double y[], double deriv[], 
			     double z1, Reionization *rz1, Astrophysics *a1, Cosmology *c1,int flag);
void clumpingBubIntegrand(double m, double y[], double deriv[]);


#endif
