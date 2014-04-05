// Astronomy.h
//
// Collate various useful astronomy conversions and definitions
//
// Want to use this for things like Magnitude conversions and so on
//
//

#ifndef ASTRONOMY_H
#define ASTRONOMY_H

#include "dcosmology.h"

//list of optical bands
enum MAGBANDS {U,B,V,R,I,J,H,K,L,M,N};
const double deltawave[11]={660.,940.,850.,1600.,1490.,2130.,3070.,3900.,4720.,4600.,0.0};
const double effwave[11]={3670.,4360.,5450.,6380.,7970.,12200.,16300.,21900.,34500.,47500.,106000.};
const double bandzero[11]={-20.94,-20.45,-21.12,-21.61,-22.27,-23.80,-24.80,-26.00,-27.87,0.0,0.0};

class Astronomy {

 public:
  Astronomy();
  ~Astronomy();

  //magnitude conversions
  double magFromL(double L, double z=0.0, double dL=1.0e-5*MPC);
  double lumFromMag(double m, double z=0.0, double dL=1.0e-5*MPC);

  double absMagFromL(double L, double z, double dL);
  double lumFromAbsMag(double M, double z, double dL);
  double absMagFromM(double M, double z, double dL);

  //Observational bands and magnitude systems
  double fluxFromMagnitude(double m, MAGBANDS band);
  double lumFromFlux(double flux, double dL, double z);
  double lumNuFromMagnitude(double m, MAGBANDS band, double dL=1.0e-5, double z=0.0); //default values assume m is absolute magnitude
  double lumBandFromMagnitude(double m, MAGBANDS band, double dL=1.0e-5, double z=0.0);

 private:

};

#endif
