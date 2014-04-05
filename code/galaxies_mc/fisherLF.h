//fisherLF.h
// information for calculating the fisher matrix for a galaxy survey
// FisherLF is a subclass of Fisher
//
// Calculates Fisher matrix for observations at a single specified redshift
// 
// Extra parameters: bias  - one per redshift : cannot connect between 
//						redshift bins
//

#ifndef FISHERLF_H
#define FISHERLF_H

#include "spline.h"
#include "observation.h"
#include "fisher.h"
#include "galaxies.h"

using namespace std;

// FisherLF class declaration

class FisherLF : public Fisher
{

 public:
FisherLF();
~FisherLF();

//initialisation
void initExperiment(string name, double zmin, double zmax, double anglex1, double angley1, double maglim1, int nfield1, int nbin1);

 //member functions

 //set reionization parameters
 void setLFParam(double phi0, double M0, double alpha);

 //spline setting
 double setGlobalFiles();
 double setGlobalPair(int iflag);

 //fisherMatrix
 void fisherMatrix();
 double getVariance(double mi, double mj);
 double getSampleVariance(double mi, double mj, CosmParam *c=NULL);
 double deriveS(double mi, double mj,CosmParam c_p,CosmParam c_m,string tag,int iflag);
 double deriveNbin(double m,CosmParam c_p,CosmParam c_m,string tag,int iflag);

 //sky temperature
 double Nbin(double m, CosmParam *c, Spline *spline=NULL);

//window
 double getWindow(CosmParam *c);

//Tinker bias
 double biasTinker(double M, double z, CosmParam *c);

//abundance matching
 double biasM(double M, double dM, double z, CosmParam *c);
 double getHaloMassFromMag(double M, double z, CosmParam *c);

protected:

 double nbin;

 double volume;
 double magbin;
 double min_mag;
 double max_mag;
 double anglex, angley;
 double min_z, max_z;
 int nfield;

 Spline splineFid;
 Spline splinePow1p;
 Spline splinePow1m;

 int model_flag;  //use 21 cm signal model
 int reset_flag;
};

//mag to mass conversions
double dummyGetMass(double mass);
double setDummyGetMass(double mass, Cosmology *c1, int iflag);
#endif
