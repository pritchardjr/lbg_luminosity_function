//fisherCMB.h
// information for calculating the fisher matrix for a CMB experiment
// FisherCMB is a subclass of Fisher

#ifndef FISHERCMB_H
#define FISHERCMB_H

#include "spline.h"
#include "observation.h"
#include "fisher.h"
using namespace std;


//const int NPRECISION(14); //sigfig for file io
const int nparamCMBX(3);

// FisherCMB class declaration

class FisherCMB : public Fisher
{

 public:
FisherCMB();
~FisherCMB();

//initialisation
void initExperiment(string name1, double fsky1, int xmax1, int lmax1, double nu1, double beam1, double sigmat1, double sigmap1);
void addChannelExperiment(double nu1, double beam, double sigmat, double sigmap);
int replaceChannel(int nreplace, double nu1, double beam1, double sigmat1, double sigmap1);
void changePolarization(int xmax1);

//member functions
void setLMax(int lmax1);

//CMB Fisher Matrix
//void derivCMB(int lmax, double p_step,CosmParam *cm, CosmParam *cp,string file);
void fisherMatrixCMB();
double formFisherCMB(int luse,double *Cl,double *dCli,double *dClj);
double covCMB(int l,int x,int y,double *Cl);

 protected:

  //CMB observation parameters
//  string exp_name;
  double fsky;
  int nchannel;
  int xmax;
  int lmax;
  double *beam;
  double *sigmat;
  double *sigmap;

};

#endif



