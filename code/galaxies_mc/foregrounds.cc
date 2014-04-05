
#include <math.h>
#include <iostream>
#include <fstream>
#include "observation.h"
#include "dcosmology.h"
#include "dnumrecipes.h"
#include "astrophysics.h"
#include "twentyonecm.h"
#include "foregrounds.h"


using namespace std;

////////////////////////////////////////////////////////////////////////
//Constructor/Destructor
/////////////////////////////////////////////////////////////////////
Foregrounds::Foregrounds(Cosmology *c1)
{
  // cout <<"21cm constructor has been called." <<endl;
	c=c1;

  //cout << "21cm Constructor called successfully." <<endl; 
}

Foregrounds::~Foregrounds()
{
  // cout <<"Atomic destructor has been called." <<endl;
}

//////////////////////////////////////////////////////////////////////////
// Member Functions
/////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// Foreground models
/////////////////////////////////////////////////////////////////////////

//Galactic synchrotron
double Foregrounds::clGalacticSynch(int l, double nu1, double nu2)
{
	double beta(2.5);
	int l0(5);
	double cl;

	cl=pow((double)(l/l0),2.0-beta);
	cl*=TSynch(nu1)*TSynch(nu2);
	cl*=2.0*PI/(double)(l*l);

	return cl;
}

// Pivot temperature for galactic synchrotron emission 
// Units:  temp  K
//	   nu    Hz
double Foregrounds::TSynch(double nu)
{
	double temp,spec;
	double Asynch(25.0);
	double nu0(150.0e6);
	double alpha(2.55);
	double dalpha(0.1);

	temp=Asynch;
	spec=-alpha-dalpha*log(nu/nu0);
	temp*=pow(nu/nu0,spec);
	
	return temp;
}


//////////////////////////////////////////////////////////////////////////
// Foreground removal
/////////////////////////////////////////////////////////////////////////
// Follows McQuinn et al. (2006)
//

double Foregrounds::legendreBasis(double nu, int order)
{
	double x,Pl;

	x=2.0*(nu-nubandmin)/bandwidth-1.0;
	Pl=gsl_sf_legendre_Pl(order,x);
	Pl*=sqrt((2.0*(double)(order)+1.0)/2.0);
	return Pl;
}
