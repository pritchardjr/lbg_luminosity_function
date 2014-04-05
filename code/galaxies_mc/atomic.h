// atomic.h
// Contains prototypes and constants needed for calculating atomic
// levels for spectroscopy calculations

#ifndef ATOMIC_H
#define ATOMIC_H

// Hated cgs units used for ease of connection to astronomy texts
const double BOHRA0(0.52918e-8);     // Bohr radius in cm 
//const double ELECTRONCHARGE(4.8033e-10);   //Electron charge in esu units
const double ELECTRONMASS(9.10938e-28);              //Electron mass in cgs
const double PLANCKCONSTANT(6.62607e-27);     // h in cgs
const double RYDBERG_CGS(2.1798719e-11);    // Rydberg in cgs
const double RYDBERG(13.60569172);          // Rydberg in Electron volts
const double SPEEDOFLIGHT_CGS(2.99792458e10);  //Speed of light in cgs
const double ELECTRONCHARGE_CGS(4.8033e-10);   //Electron charge in esu units

static double **pstore;

class Atomics {

 public:
  Atomics(int rNMAX_in, double **pastore, double **pastore_thick);
  ~Atomics();

  // Utility functions

  // Hydrogen Transition Probability function
  double transitionS(double n, double l, double j, double np, double lp,
		     double jp);
  double transitionA(double n, double l, double j, double np, double lp,
		     double jp);
  double transitionATotal(double n, double l, double j);
  double hydrogenFreq(double n, double np);
  double hydrogenWave(double n, double np);
  double hydrogenEnergy(double n, double np);
  double hydrogenMatrix(double n, double l, double np, double lp);
  double polyR(double n, double l, double np, double lp);
  double transitionP(double n, double l, double j, double np, double lp,
		     double jp);
  double transitionF(double n, double l, double j, double np, double lp,
		     double jp);

  double recycleLyman(double n, double l, double j, int nmax);
  double initRecycleLyman(double **pastore);
  double initRecycleLymanThick(double **pastore_thick);

  //  Math functions for Hypergeometric Functions
  double intPochhammer(double a, double k);
  double term2F1(double a, double b, double c, double x);
  double factorial(double n);

  // Data retreival functions
  double getRecycleLyman(double n, double l, double j, double **pastore);
  double getRecycleLymanThick(double n, double l, double j, double **pastore_thick);


 protected:
  double ryd;
  int rNMAX;
  int thick_flag;



};


#endif
