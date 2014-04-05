// observation.h
// Specifies details of 21 cm observations and their sensitivity
// Also handle foregrounds here, so basically the practicalities.

#ifndef OBSERVATION_H
#define OBSERVATION_H

#include "dcosmology.h"
#include "spline.h"

const int MCQUINN(0);
const int NPRECISION(14); //sigfig for file io

class Observation {

 public:
  Observation(Cosmology *c1);
  ~Observation();

//member functions
void setCosmology(Cosmology *c1);
double getDMin();
double getDMax();
double getAlpha();
double getNAnt();
double getFreqRes();
double getNInnerRing();
double getRInnerRing();

// fisher integration limits
double getThetaMin(double z, double k);
double getThetaMax(double z, double k);
double getSampleConstantD(double z);
double getThermalConstantE(double z);

//Observational performance
double fieldOfView(double nu);
double angularWidth(double nu);
double angularResolution(double nu);

double surveyDepth(double z);
double surveyWidth(double z);
double surveyVolume(double z);
double surveyDistance(double z);

double getAntennaeArea();
double effectiveAntennaeArea(double nu);

double foregroundCutoff(double z);

// Antennae distribution
double antennaeDensity(double r);
double findDensityBreak();
double antennaeDensityFull(double r);
double findDensityBreakFull();
void optimiseAntennae(double z);

//Visibility distibution
double xcolAntennae(double r);
double baselineDensityUV(double u, double z);
void initXColSpline();
double getXCol(double r);
double numberVisibilities();
void forceNormVisibilities();


// Antennae response function
double antennaAngularResponse(double theta, double nu);

//Observation initialisation
	void initObservation(double Nant1, double Atot1, double Dmin1, double Dmax1, double B1, double Tint1, double eta1, double alpha1, double Nfield1, double antennaeProp[]);
	void initObservationFromFile(string file);
	void saveObservationToFile(string file);

//Power spectrum sensitivities
	double tempSky(double nu);
	double powerSpectrumError(double z, double k, double ps);
	double powerSpectrumErrorMu(double z, double k, double ps0, double ps2, double ps4, int order);
	double fisherErrorMu(double z, double k, double ps0, double ps2, double ps4, int order);

 protected:
	Cosmology *c;
// Base survey parameters
	double Nant;    //Total number of antennae
	double Atot;    //total effective collecting area at z=8
	double Dmin;    //min baseline=size of antenna
	double Dmax;    //max baseline
	double B;       //Bandwidth
	double Tint;    //integration time
	double eta;     //logarithmic step size dk=eta*k
	double alpha;   //antennae distribution index
	double freqRes; //frequency resolution
	double Nfield;  //number of fields observed

	Spline visDens; //visibility density distribution
	double xColRenorm;

	int antennaeDistFlag;
	double rbreakDist;
	double rInnerRing;
	double rOuterRing1;
	double rOuterRing2;
	double nInnerRing;
	double nOuterRing1;
	double nOuterRing2;

	double ZOPTIMISE;
	int foregroundFlag;
	int optimiseFlag;
  	int fixed_area_flag;
};

// Antenna distribution
double setDensityBreak(double r, Observation *Obs1, int iflag);
double dummyDensityBreak(double r);
double setDensityBreakFull(double r, Observation *Obs1, int iflag);
double dummyDensityBreakFull(double r);

//visibility distribution
double dummyRXCol(double r);
double setRXCol(double r, Observation *Obs1, int iflag);
double dummyPhiXCol(double phi);
double setXColAntennae(double rp, double phi, Observation *Obs, int iflag);
double setNumberVisibilitiesInt(double r, Observation *Obs1, int iflag);
double dummyNumberVisibilitiesInt(double r);


//power spectrum error kernel
double setPowerSpectrumErrorKernel(double theta, double k1, double z1, double x1, double D1, double E1, double ps1, Observation *Obs1, int iflag);
double dummyPowerSpectrumErrorKernel(double theta);
//power spectrum error kernel - mu
double setFisherErrorMuKernel(double theta, double k1, double z1, double x1, double D1, double E1, double ps01, double ps21, double ps41, int order1, Observation *Obs1, int iflag);
double dummyFisherErrorMuKernel(double theta);

//matrix manipulation
void invertMatrix(double **fisher,int nparam, double **ifisher);
void invertMatrixZO(double **fisher,int nparam, double **ifisher);
void gaussj(double **a, int n, double **b, int m);
void writeMatrix(double **M,int nparam, char *file);
void writeMatrixS(double **M,int nparam, string file);
double detMatrix(double **M, int nparam);
void multMatrix(double **A, double **B, double **C, int N);
void invertMatrix2(double **fisher,int nparam, double **ifisher);
void invertMatrixSVD(double **fisher,int nparam, double **ifisher);
#endif
