// Miscellaneous programs for my radio project
//
//


#include <math.h>
#include <iostream>
#include <fstream>
#include "astrophysics.h"
#include "dcosmology.h"
#include "dnumrecipes.h"
#include "spline.h"
#include "radio.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////
// Constructor/Destructor 
///////////////////////////////////////////////////////////////////////////
Radio::Radio(Cosmology *c1, Astrophysics *a1)
{
  c=c1;
  a=a1;
}

Radio::~Radio()
{

}





////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


//assume flat Universe
double dndzdo(double z, Cosmology *c)
{
  double dn,zz;
  double r,cH;

  zz=z+1.0;
  cH=SPEEDOFLIGHT_CGS/H0_CMSMPC/c->getH();
  r=c->coordDistance(z)*cH;
  dn=c->nCollObject(z,coolMass(c,z));    //comoving number density
  dn*=cH/c->hubbleZ(z);
  dn*=r*r;

  return dn;
}

double dndo(double z1, double z2, Cosmology *c)
{
  double dn;
  double tol(1.0e-4);

  setDnDo(0.0,c,1);
  dn=qromb(getDnDo,z1,z2,tol);

  return dn;
}

double setDnDo(double z, Cosmology *c1, int iflag)
{
  static Cosmology *c;
  double dn;

  if(iflag==1){
    c=c1;
    return 0.0;
  }

  dn=dndzdo(z,c);

  return dn;
}

double getDnDo(double z)
{
  return setDnDo(z,NULL,0);
}
