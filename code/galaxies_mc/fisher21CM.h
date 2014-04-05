//fisher21CM.h
// information for calculating the fisher matrix for a 21 cm survey
// Fisher21CM is a subclass of Fisher
//
// Calculates Fisher matrix for observations at a single specified redshift
// 
// Extra parameters: bias  - one per redshift : cannot connect between 
//						redshift bins
//

#ifndef FISHER21CM_H
#define FISHER21CM_H

#include "spline.h"
#include "observation.h"
#include "fisher.h"
#include "dcosmology.h"
using namespace std;


// Fisher21CM class declaration

class Fisher21CM : public Fisher
{

 public:
Fisher21CM();
~Fisher21CM();

//initialisation routines
void initObservation(string name1, double z1, double kmin1, double kmax1, Observation *Obs1, int regime_flag);

// member functions
string getSurveyName();
double getSurveyZ();
double getKMin();
double getKMax();
void setKMax(double kmax1);

//bias parameters
double getBias(CosmParam *c);

//set useful parameters
void setSurveyConstants();
void setNoThermalNoise();

//set galaxy files
double setGalaxyFiles();
double setGalaxyPair(int iflag);

//calculate fisher matrix
void fisherMatrix21CM();

// calculate experimental uncertainty
double variance21CM(double k, double theta);
double effectiveVolume(double k, double theta);

//calculate derivatives of power spectrum
double derivePower21CM(double k, double theta, CosmParam c_p, CosmParam c_m, string tag, int iflag);
void setDerivSplines(string tag1, string tag2, int iflag);

//power spectrum code
double getPowerSpectrumObs(double kref, double muref,CosmParam *c, Spline *spline);
//double powerSpectrum(double k, CosmParam *c, Spline *spline);

//errors on k power spectrum
double powerSpectrumError(double k);
void outPowerErrors();
double averagedEffectiveVolume(double k);

protected:

  Observation *SurveyObs;

  //survey parameters
  double survey_z;
  double kmin;
  double kmax;

  double sampleConstantD;
  double thermalConstantE;
  double xSurveyVolume;  //comoving distance to survey volume

  int post_reion_flag;

  //splines
  Spline splineFid;
  Spline splinePow1p;
  Spline splinePow1m;
  Spline splinePow2p;
  Spline splinePow2m;

};

//fisher calculation
void setFisherKernelK21CM(double lk, double y[], double deriv[],double z1, Observation *Obs1, int iflag);
void getFisherKernelK21CM(double lk, double y[], double deriv[]);
void getFisherKernel21CM(double theta, double y[], double deriv[]);
void setFisherKernel21CM(double theta, double y[], double deriv[], string file_fid, Fisher21CM *survey1,string file_dlPi1, string file_dlPj1,int iflag);

double setKErrorKernel21CM(double k1, double theta, Fisher21CM *survey1, int iflag);
double getKErrorKernel21CM(double theta);

double setKVolumeKernel21CM(double k1, double theta, Fisher21CM *survey1, int iflag);
double getKVolumeKernel21CM(double theta);

#endif