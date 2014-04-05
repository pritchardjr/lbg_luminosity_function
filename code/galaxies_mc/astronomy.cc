/*  astronomy.cc
 *
 *  Stuff for calculating astronomy related conversions and quantities
 *
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include "astronomy.h"
#include "dnumrecipes.h"  

using namespace std;

//References:

/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

Astronomy::Astronomy()
{

}

Astronomy::~Astronomy()
{
  // cout <<"Reionization destructor has been called." <<endl;

}

////////////////////////////////////////////////////////////////////
// AB magnitude conversions
////////////////////////////////////////////////////////////////////

// Calculates the apparent AB magnitude given the luminosity
//Units: L  erg s^-1 Hz^-1
//       flux erg s^-2 cm^-2 Hz^-1
//       dL cm
double Astronomy::magFromL(double L, double z, double dL)
{
  double flux;
  double m;

  //calculate flux assuming luminosity in a freq. interval
  flux=L*(1.0+z)/(4.0*PI*dL*dL);
  m= -2.5*log10(flux)-48.60;

  return m;
}

// Calculates the luminosity given the apparent AB magnitude
//Units: L  erg s^-1 Hz^-1
double Astronomy::lumFromMag(double m, double z,double dL)
{
  double flux;
  double L;
  L=exp(-0.4*(m+48.60)*log(10.0));
  L*=4.0*PI*dL*dL/(1.0+z);

  return L;
}

// Calculates the absolute AB magnitude given the luminosity
//Units: L  erg s^-1 Hz^-1
double Astronomy::absMagFromL(double L, double z, double dL)
{
  double M;
  M=magFromL(L,z,dL)+2.5*log10(1+z)-5.0*log10(dL)+5.0;
  return M;
}

// Calculates the absolute AB magnitude given the relative mag
//Units: L  erg s^-1 Hz^-1
//Note: expects dL in Parsecs
double Astronomy::absMagFromM(double m, double z, double dL)
{
  double M;
  M=m+2.5*log10(1+z)-5.0*log10(dL)+5.0;
  return M;
}

// Calculates the luminosity given the absolute AB magnitude
//Units: L  erg s^-1 Hz^-1
//
// Note: expects dL in Parsecs
double Astronomy::lumFromAbsMag(double M, double z, double dL)
{
  double m, L;

  //convert M to apparent magnitude
  m=M-2.5*log10(1+z)+5.0*log10(dL)-5.0;
//NOTE think I may be confused over luminosity and 
// luminosity/Hz in this expressions

  //use conversion from apparent mag to L
  L=lumFromMag(m,z,dL);

  cout<<"L1="<<L<<endl;
  //simpler conversion
  //M= -2.5*log10(L)+5.0+2.5*log10(4.0*PI)-48.60;
  
  L= -M-48.60+5.0+2.5*log10(4.0*PI);
  L=pow(10.0,L/2.5);
  cout<<"L2="<<L<<endl;

  return L;
}

////////////////////////////////////////////////////////////////////////
// Convert different magnitude systems to a common standard
////////////////////////////////////////////////////////////////////////
//enum MAGBANDS {U,B,V,R,I,J,H,K,L,M,N};

// Given apparent magnitude m for an observed band calculate the corresponding
// flux that the detector sees.
//
// Units: m is apparent magnitude
//        F is the flux in erg/s/cm^2/Hz
//
double Astronomy::fluxFromMagnitude(double m, MAGBANDS band)
{
   double F;
   double zeropnt;
   double lambda;

   //Johnson-Morgan bands based on A0V stars
   //fom webpage
   //double effwave[11]={3600.,4400.,5500.,6600.,8000.,12500.,16000.,21800.,34500.,47500.,106000.};
   //double bandzero[11]={-20.94,-20.45,-21.12,-21.61,-22.27,-23.80,-24.80,-26.00,-27.87,0.0,0.0};

   //from Schneider - more accurate band details with the effective wavelength 
   //and band width in Angstroms
   //double deltawave[11]={660.,940.,850.,1600.,1490.,2130.,3070.,3900.,4720.,4600.,0.0};
   //double effwave[11]={3670.,4360.,5450.,6380.,7970.,12200.,16300.,21900.,34500.,47500.,106000.};

   //magnitude zeropnt is for a source having flux 1 erg/s/cm^2/Angstrom
   //this is not the same as erg/s/cm^2/Hz !!!

   //difference to convert between these units corresponds to a factor of
   // Fnu=Flam*lambda^2/c
   // (Angstrom/Hz)=(1e-10/c) multiplying the flux or a difference in 
   // magnitude of delta m= 46.1928

   //correctly calculate flux as Flam in erg/s/cm^2/Angstrom
   zeropnt=bandzero[band];
   lambda=effwave[band];
   F=pow(10.0,-0.4*(m-zeropnt));

   //cout<<F<<endl;
   //Convert to Fnu in erg/s/cm^2/Hz
   //this is a little dubious, since really depends upon 
   //details of ther eference objects and filters to convert from nu to lam.
   F*=pow(lambda,2.0)*1.0e-8/SPEEDOFLIGHT_CGS;

   return F;
}

////////////////////////////////////////////////////////////////////////
// Given flux and distance calculate the source luminosity
////////////////////////////////////////////////////////////////////////

// Convert flux to luminosity
// Units: flux : erg/s/cm^2/Hz
//        dL  : Mpc
//        Lnu : erg/s/Hz
double Astronomy::lumFromFlux(double flux, double dL, double z)
{
   double Lnu, area;

   Lnu=flux;
   area=4.0*PI*dL*dL/(1.0+z)*MPC*MPC;
   Lnu*=area;

   return Lnu;
}

// Convert magnitude into band luminosity
// Units: 
//        m : apparent magnitude
//        flux : erg/s/cm^2/Hz
//        dL  : Mpc
//        Lnu : erg/s/Hz
//
double Astronomy::lumNuFromMagnitude(double m, MAGBANDS band, double dL, double z)
{
   double Fnu, Lnu;

   Fnu=fluxFromMagnitude(m,band);
   Lnu=lumFromFlux(Fnu,dL,z);
   return Lnu;
}



// Convert magnitude into band luminosity
// Here multiply by the filter width to get the total luminosity inferred
// from all flux detected in the band
// Units: 
//        m : apparent magnitude
//        flux : erg/s/cm^2/Hz
//        dL  : Mpc
//        Lband : erg/s
//
double Astronomy::lumBandFromMagnitude(double m, MAGBANDS band, double dL, double z)
{
   double Fnu, Lnu, Lband;
   double dnu, lambda;

   Fnu=fluxFromMagnitude(m,band);
   Lnu=lumFromFlux(Fnu,dL,z);
   lambda=effwave[band]*1.0e-8;
   dnu=deltawave[band]*1.0e-8*SPEEDOFLIGHT_CGS/lambda/lambda;
   Lband=Lnu*dnu;
   return Lband;
}
