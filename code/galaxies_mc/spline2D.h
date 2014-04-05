// spline2D.h
//
// Class for creating 2 dimensional splines of a table of data
// different class instances for different data sets saving 
// redundancy of code when you need to spline multiple quantities
//
// matrix is m by n with m values of x1 and n values of x2
//
// Selfcontained piece of code.  Could use some of the code in Spline class
// to make more compact, but not really worth the trouble at present.

#ifndef SPLINE2D_H
#define SPLINE2D_H

#include "dnumrecipes.h"
#include "spline.h"
using namespace std;

class Spline2D {

 public:
  Spline2D();
  ~Spline2D();

  //member function
  int getN();
  int getM();
  double getX1Max();
  double getX1Min();
  double getX2Max();
  double getX2Min();
  double getX1element(int i);
  double getX2element(int i);
  int checkInit();

  //
  void setSplineSP(int m1, int n1, double x1[], double x2[], double **y);
  void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
  void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y);


  void splineSP(double x[], double y[], int n, double yp1, double ypn,
		double y2[]);
  void splintSP(double xa[], double ya[], double y2a[], int n, double x,
	      double *y);
  double returnValue(double x1, double x2);
  void hunt(double xx[], int n, double x, int *jlo);

  void loadFileSpline(string file);
  void cleanSpline();

 protected:
  double *x1as,*x2as, **yas, **y2as;
  double x1min,x1max,x2min,x2max;
  int n, m;
  int init_flag;

};


#endif
