// Specifies details of 21 cm foregrounds


#ifndef FOREGROUNDS_H
#define FOREGROUNDS_H

#include "dcosmology.h"
#include "spline.h"
#include "observation.h"
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_legendre.h>

class Foregrounds {

public:
	Foregrounds(Cosmology *c1);
	~Foregrounds();

double clGalacticSynch(int l, double nu1, double nu2);
double TSynch(double nu);

//foreground removal
double legendreBasis(double nu, int order);

protected:
	Cosmology *c;
	double nubandmin;     //lowest frequency of band
	double nubandmax;     //highest frequency of band
	double bandwidth;     //Observation bandwidth
	double deltanu;       //frequency resolution of bins

	int npolyfit;  	      //order of polynomial to fit to data
	int nbands;           //number of frequency bands

};


#endif