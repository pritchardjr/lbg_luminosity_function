/* haloDensity.cc
 * 
 * Program for estimating the density power spectrum with the halo
 * model.  The formalism is described in Cooray & Sheth (2002).
 */
#include <iostream>
#include "dnumrecipes.h"
#include "math.h"
#include "dcosmology.h"
#include "haloDensity.h"

using namespace std;

/********************************************************************
 ********************* Power Spectrum of Density ********************
 ********************************************************************/

/* Normalized density profile of an NFW halo. 
 *   r = distance from center
 *   c = halo concentration parameter
 *   rvir = virial radius */
double urhoNFW(double r, double c, double rvir)
{
  double fc,cx;
  double ans;

  fc = (log(1.0+c)-c/(1.0+c));
  cx = c*r/rvir;
  ans = 1.0/(4.0*PI*fc)*(c*c*c)/(rvir/rvir/rvir);
  ans *= 1.0/(cx*(1.0+cx)*(1.0+cx));
  return ans;
}

/* Fourier transform of normalized density profile of an NFW halo. 
 *   k = wavenumber (in comoving Mpc^-1) 
 *   c = halo concentration 
 *   rs = break radius (in comoving Mpc) */
double uNFW(double k, double c, double rs)
{
  double kr;
  double krc;
  double sic,si;
  double cic,ci;
  double fc;
  double ans;
  double sidiff,cidiff;

  kr = (k*rs);
  krc = (kr*(1.0+c));

  fc = (log(1.0+c)-c/(1.0+c));
  cisi(kr,&ci,&si);
  cisi(krc,&cic,&sic);
  sidiff = sic - si;
  cidiff = cic - ci;
  //  cout << "sidiff=" << sidiff << " cidiff=" << cidiff << endl;
  ans = sin(kr)*sidiff/fc - sin(kr*c)/krc/fc + cos(kr)*cidiff/fc;
  return (double)(ans);
}

/* Returns concentration of a halo, using the Bullock et al. fit 
 *   m = halo mass (Msun) */
double concentration(double m, double z)
{
  double mstar(1.0e9);
  double p(-0.1);

  return (25.0/(1.0+z)*pow(m/mstar,p));
}

/* Returns break radius of a halo in comoving Mpc.
 *   m = halo mass (Msun)
 *   conc = halo concentration */
double findRs(double m, double conc, double z, Cosmology *c)
{
  double rv;

  rv = 1.0e-3*rvir(c,m,z)*(1.0+z);
  return (rv/conc);
}

/* Compute the one-halo contribution to the density power spectrum by 
 * integrating over the halo mass function.  Uses setPS1hddIntegrand
 * as the integrand.
 *   k = wavenumber (comoving Mpc^-1) */
double PS1hdd(double k, double z, Cosmology *c)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double ps;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setPS1hddIntegrand(0.0,ans,ans,k,z,c,1);
  mMin = filterMass(c,z);
  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
	   ps1hddIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    mMin = mMax;
  }
  ps = ans[1]/(CRITDENMSOLMPC*CRITDENMSOLMPC*c->getOm0hh()*c->getOm0hh());
  free_dvector(ans,1,1);
  return ps;
}

/* Integrand for computing the 1-halo term in the power spectrum.
 *   m = halo mass (Msun) [integration variable]
 *   y,deriv = parameters used in the integration
 *   k1 = wavenumber (comoving Mpc^-1)
 */
void setPS1hddIntegrand(double m, double y[], double deriv[], double k1, 
			double z1, Cosmology *c1, int flag)
{
  static double k,z;
  static Cosmology *c;
  double con,rs;
  double u;
  double numberDensity;

  if (flag == 1) {
    k = k1;
    z = z1;
    c = c1;
    return;
  }

  con = concentration(m,z);
  rs = findRs(m,con,z,c);
  u = uNFW(k,con,rs);
  //  numberDensity = c->dndlM(z,m)/m;
  numberDensity = nm(m,z,c)/m;
  deriv[1] = numberDensity*m*m*u*u;
  /*
  cout << m << "\t" << vol << "\t" << chi << "\t" << numberDensity 
       << "\t" << deriv[1] << endl;
  */
}

/* Dummy function for integrating the 1-halo term */
void ps1hddIntegrand(double m, double y[], double deriv[])
{
  setPS1hddIntegrand(m,y,deriv,0.0,0.0,NULL,0);
}

/* Compute the 2-halo contribution to the density power spectrum,
 * assuming that you already know the necessary integral and Plin. 
 *   integral = the 2-halo prefactor, for the part where we explicitly
 *              integrate over the mass function
 *   renorm = At high redshifts, most of the mass is in tiny halos which 
 *            don't converge well.  The factor renorm, to be computed with 
 *            RenormInt, adds the contributions of these halos back in.
 *   Plin = linear PS at this wavenumber. */
double PS2hdd(double integral, double renorm, double Plin)
{
  return ((integral+renorm)*(integral+renorm)*Plin);
}

/* Compute the 2-halo contribution to the density power spectrum, assuming that
 * you haven't precomputed the necessary integral. 
 *   k = wavenumber (comoving Mpc^-1)
 *   renorm = see PS2hdd
 */
double PS2hddFull(double k, double z, double renorm, Cosmology *c)
{
  double ps,integral;

  integral = PS2hddInt(k,z,c) + renorm;
  ps = integral*integral;
  ps *= Plinear(k,z,c);
  return ps;
}

/* Perform the necessary integral for the 2-halo contribution to the
 * density power spectrum. The integrand is setPS2hddIntegrand.  We
 * will only compute the integral explicitly for halos above the
 * filter mass; this is fine at low z, but at high z a large fraction
 * of the mass is contained in smaller halos.  These are accounted for
 * in renormInt below.
 *   k = wavenumber (comoving Mpc^-1) */
double PS2hddInt(double k, double z, Cosmology *c)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double ps;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setPS2hddIntegrand(0.0,ans,ans,k,z,c,1);
  mMin = filterMass(c,z);
  //  mMin = 1.0;
  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
	   ps2hddIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    mMin = mMax;
  }
  ps = ans[1]/(CRITDENMSOLMPC*c->getOm0hh());
  free_dvector(ans,1,1);
  return ps;
}

/* Integrand for prefactor to 2-halo term.
 *   m = halo mass (Msun) [integration variable]
 *   y,deriv = parameters used in the integration
 *   k1 = wavenumber (comoving Mpc^-1)
 */
void setPS2hddIntegrand(double m, double y[], double deriv[], double k1, 
			double z1, Cosmology *c1, int flag)
{
  static double k,z;
  static Cosmology *c;
  double rs,con;
  double u,bias,numberDensity;

  if (flag == 1) {
    k = k1;
    z = z1;
    c = c1;
    return;
  }

  con = concentration(m,z);
  rs = findRs(m,con,z,c);
  u = uNFW(k,con,rs);
  //  bias = c->biasPS(z,m);
  bias = biasm(m,z,c);
  //  numberDensity = c->dndlM(z,m)/m;
  numberDensity = nm(m,z,c)/m;
  deriv[1] = numberDensity*m*u*bias;
  /*
  cout << m << "\t" << con << "\t" << u << "\t" << bias << "\t" 
       << m*numberDensity << "\t" << deriv[1] << endl;
  */
}

/* Dummy function for integrating the 2-halo term */
void ps2hddIntegrand(double m, double y[], double deriv[])
{
  setPS2hddIntegrand(m,y,deriv,0.0,0.0,NULL,0);
}

/* Perform the necessary integral to add in the contribution from small
 * halos to the 2-point density power spectrum result. For these halos, 
 * we are at k much smaller than that which corresponds to the break radius,
 * so u->1 and we neglect that.  We therefore just want 
 *           \int_0^M_fil dm m n(m) b(m)
 * Computing this is difficult, because the sigma(m) gets difficult to 
 * evaluate at small m.  We therefore note that (Cooray & Sheth, eq. 71)
 *           \int_0^\infty dm m n(m) b(m) = 1 
 * and so compute the high-mass end of the integral and subtract. */
double RenormInt(double z, Cosmology *c)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double ps;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setRenormIntegrand(0.0,ans,ans,z,c,1);
  mMin = filterMass(c,z);
  //  mMin = 1.0;
  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
	   renormIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    mMin = mMax;
  }
  ps = 1.0 - ans[1]/(CRITDENMSOLMPC*c->getOm0hh());
  free_dvector(ans,1,1);
  return ps;
}

/* Integrand for computing the small-halo contribution to the 2-halo term.
 *   m = halo mass (Msun) [integration variable]
 *   y,deriv = parameters used in the integration
 */
void setRenormIntegrand(double m, double y[], double deriv[], double z1,
			Cosmology *c1, int flag)
{
  static double z;
  static Cosmology *c;
  double bias,numberDensity;

  if (flag == 1) {
    z = z1;
    c = c1;
    return;
  }

  //  bias = c->biasPS(z,m);
  //  numberDensity = c->dndlM(z,m)/m;
  bias = biasm(m,z,c);
  numberDensity = nm(m,z,c)/m;
  deriv[1] = numberDensity*m*bias;
  //  cout << m << "\t" << bias << "\t" 
  //       << m*numberDensity << "\t" << deriv[1] << endl;
}

/* Dummy variable for the small-halo contribution */
void renormIntegrand(double m, double y[], double deriv[])
{
  setRenormIntegrand(m,y,deriv,0.0,NULL,0);
}

/********************************************************************
 ***************** Correlation Function of Density ******************
 ********************************************************************/

/* Construct an interpolation table for the density correlation function.
 *   rad = comoving distance (Mpc)
 * If flag=1, sets up the interpolation table (if the appropriate file 
 * exists, reads it in; otherwise computes the table and writes the file).  
 * If flag=2, clears table from memory.  If flag=0, interpolates.
 * Uses spline interpolation from Numerical Recipes.
 * To set up the \xi(r) interpolation table, it first computes a P(k)
 * table and then takes the (3D, isotropic) Fourier Transform.
 */
double xidd(double rad, double z, Cosmology *c, int flag)
{
  double natural(1.0e30);
  double *ans,oldAns;
  int goodSteps,badSteps;
  double step;
  static double rMin,rMax;
  //  const int nrad(220);
  int nrad(200);
  static double *r;
  static double *xivals,*xisp;
  double deltar,lgrmin,lgrmax;
  double xi;
  double kMin,kMax;
  int i;
  char xiddFile[100];
  FILE *infile; 

  //nrad=400;

  if (flag == 1) {

    r = dvector(1,nrad);
    xivals = dvector(1,nrad);
    
    /* Initialize the power spectrum interpolation table */
    PSddTable(0.0,z,c,1);

    /* Construct the filename */
    sprintf(xiddFile,"/Users/jon/Documents/projects/include/xidd/xidd_z%05.2lf.dat",
	    z);

    /* Check if the file exists */
    if(!(infile=fopen(xiddFile, "r"))) {
      /* We need to make the data (and print it out) */
      infile=fopen(xiddFile, "w");
      cerr << "\nGenerating file..." << xiddFile << endl;
     
      /* Set up the wavenumber dvector */
      rMin = 1.0e-4;
      rMax = 2.5e3;
      //  rMin = 1.0e-4;
      //rMax = 8.0e3;
      //      rMax = 3.0e3;
      //rMax=1.0e4;

      lgrmin = log(rMin);
      lgrmax = log(rMax);
      deltar = (lgrmax-lgrmin)/(double)(nrad);

      for (i=1;i<=nrad;i++) {
	r[i] = exp(lgrmin + (double)(i-1)*deltar);
      }

      ans = dvector(1,1);
      for (i=1;i<=nrad;i++) {
	ans[1] = 0.0;
	setCorrInt(0.0,ans,ans,r[i],z,c,1);
	kMin = 1.0e-6;
	while (1) {
	  kMax = kMin + 2.0*PI/r[i];
	  if (kMax > kMin*10.0)
	    kMax = kMin*10.0;
	  step = (kMax-kMin)/10.0;
	  oldAns = ans[1];
	  odeint(ans,1,kMin,kMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
		 corrInt,bsstep);
	  //	  cout << "kMin=" << kMin << "\tkMax=" << kMax << "\tans="
	  //	       << ans[1] << endl;
	  if (fabs((ans[1]-oldAns)/ans[1]) <= 1.0e-6)
	    break;
	  kMin = kMax;
	}
	xivals[i] = ans[1];
	fprintf(infile,"%g %g\n",r[i],xivals[i]);
      }
      free_dvector(ans,1,1);

      fclose(infile);
    } else {
      cerr << "\nUsing file... " << xiddFile << endl;
      /* File exists, just read it in */
      for (i=1;i<=nrad;i++) {
	fscanf(infile,"%lf %lf\n",&r[i],&xivals[i]);
      }
      fclose(infile);

      rMin = r[1];
      rMax = r[nrad];
    }

    /* Set up interpolation table. */
    xisp = dvector(1,nrad);

    spline(r,xivals,nrad,natural,natural,xisp);

    return 0.0;

  } else if (flag == 2) {

    free_dvector(xisp,1,nrad);
    free_dvector(xivals,1,nrad);
    free_dvector(r,1,nrad);
    PSddTable(0.0,z,c,2);

    return 0.0;
  }

  if (rad >= rMin && rad <= rMax) {
    splint(r,xivals,xisp,nrad,rad,&xi);
  } else if (rad < rMin) {
    /* Plug in the minimum value in this regime (okay b/c the 
     *  volume will be so small anyway) */
    xi = xivals[1];
  } else {
    ans = dvector(1,1);
    ans[1] = 0.0;
    setCorrInt(0.0,ans,ans,rad,z,c,1);
    kMin = 1.0e-6;
    while (1) {
	kMax = kMin + 2.0*PI/rad;
	if (kMax > kMin*10.0)
	  kMax = kMin*10.0;
	step = (kMax-kMin)/10.0;;
	oldAns = ans[1];
	odeint(ans,1,kMin,kMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
	       corrInt,bsstep);
	//	cout << "kMin=" << kMin << "\tkMax=" << kMax << "\tans="
	//	     << ans[1] << endl;
	if (fabs((ans[1]-oldAns)/ans[1]) <= 1.0e-6)
	  break;
	kMin = kMax;
    }
    xi = ans[1];
    free_dvector(ans,1,1);
  }
  return xi;
}

/* Integrand to compute the correlation function, given the P(k).
 *   kp = wavenumber (comoving Mpc^-1) [integration variable]
 *   y,deriv = parameters for integration routine
 *   r1 = distance for which you are computing the correlation function
 */
void setCorrInt(double kp, double y[], double deriv[], double r1, 
		double z1, Cosmology *c1, int flag)
{
  static double rad,z;
  static Cosmology *c;
  double psdd;

  if (flag == 1) { 
    rad = r1;
    z = z1;
    c = c1;
    return;
  }

  /* Evaluate the power spectrum */
  psdd = PSddTable(kp,z,c,0);
  deriv[1] = kp*kp*psdd/(2.0*PI*PI)*sin(kp*rad)/(kp*rad);
  //  cout << kp << "\t" << psdd << "\t" << deriv[1] << endl;
  return;
}

/* Dummy integrand for the correlation function */
void corrInt(double kp, double y[], double deriv[])
{
  setCorrInt(kp,y,deriv,0.0,0.0,NULL,0);
}

/* Constructs an interpolation table for the power spectrum. (Including
 * all terms).
 * If flag=1, actually makes the table.  If flag=2, clears it from 
 * memory.  If flag=0, computes P(k)
 * Uses spline interpolation from Numerical Recipes.
 */
double PSddTable(double kp, double z, Cosmology *c, int flag)
{
  double natural(1.0e30);
  static double kMin,kMax;
  const int nk1(400);
  static double *lgk;
  static double *psdd,*psddsp;
  double k;
  double deltak,lgkmin,lgkmax;
  double lgkp,lgps,ps;
  static double renorm;
  int i;
  double pslin,ps1dd,ps2dd,dInt;

  if (flag == 1) {

    lgk = dvector(1,nk1);
    psdd = dvector(1,nk1);

    /* Set up the wavenumber dvector */
    kMin = 1.0e-6;
    kMax = 1.0e6;
    lgkmin = log(kMin);
    lgkmax = log(kMax);
    deltak = (lgkmax-lgkmin)/(double)(nk1);

    for (i=1;i<=nk1;i++) {
      lgk[i] = lgkmin + (double)(i-1)*deltak;
    }

    /* Compute the small-halo correction to the 2-halo term */
    renorm = RenormInt(z,c);

    for (i=1;i<=nk1;i++) {
      k = exp(lgk[i]);
      pslin = Plinear(k,z,c);
      ps1dd = PS1hdd(k,z,c);
      dInt = PS2hddInt(k,z,c);
      ps2dd = PS2hdd(dInt,renorm,pslin);
      psdd[i] = log((ps1dd + ps2dd + 1.0e-30));
      //      cout << k << "\t" << pslin << "\t" << ps1dd << "\t" << ps2dd 
      //      	   << "\t" << exp(psdd[i]) << endl;
    }

    /* Set up interpolation table. */
    psddsp = dvector(1,nk1);

    spline(lgk,psdd,nk1,natural,natural,psddsp);

    return 0.0;

  } else if (flag == 2) {

    free_dvector(psddsp,1,nk1);
    free_dvector(psdd,1,nk1);
    free_dvector(lgk,1,nk1);

    return 0.0;
  }

  if (kp >= kMin && kp <= kMax) {
    lgkp = log(kp);
    splint(lgk,psdd,psddsp,nk1,lgkp,&lgps);
    ps = exp(lgps);
  } else {
    //    cout << "Out of bounds in PSTable!!! kp = " << kp << "\n";
    pslin = Plinear(kp,z,c);
    ps1dd = PS1hdd(kp,z,c);
    dInt = PS2hddInt(kp,z,c);
    ps2dd = PS2hdd(dInt,renorm,pslin);
    ps = ps1dd + ps2dd;
  }
  return ps;
}

/********************************************************************
 ************************* Utility Functions ************************
 ********************************************************************/

/* Compute the linear power spectrum at redshift z.  Uses the Eisenstein
 * & Hu fit to the transfer function. 
 *   k = comoving wavenumber (Mpc^-1) */
double Plinear(double k, double z, Cosmology *c)
{
  return (c->powerSpectrum(k)*pow(c->growthFac(z),(double)(2.0)));
}


/////////////////////////////////////////////////////////////////////
// Jonathan's additions to Steve's Code
/////////////////////////////////////////////////////////////////////

/* Halo bias (linear) */
//
// The expression here is taken from Mo and White (1996), which agrees with
// Cooray and Sheth 2002.  There's a difference of the growth function in
// the denominator of the second term.  Not sure how to explain this.
//
//
double haloBias(double mass, double z, Cosmology *c)
{
  double deltasc,deltasc0,sig;
  double bias;
  double dummy;

  deltasc = c->delCrit0(z)/c->growthFac(z);
  deltasc0 = c->delCrit0(0.0);
  sig = sigm(mass,dummy,c,0);
  bias = 1.0 + (deltasc*deltasc/sig/sig - 1.0)/deltasc0;
  return bias;
}

//mass averaged bias for halos
double meanHaloBias(double z, Cosmology *c)
{
    double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setHaloBiasIntegrand(0.0,ans,ans,z,c,1);
  mMin=coolMass(c,z);  //minimum halo mass

  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];

    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
    	   haloBiasIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    mMin = mMax;
  }
  bias = ans[1]/(CRITDENMSOLMPC*c->getOm0hh());
  bias /= c->fColl(z);  //normalise halo mass function appropriately
  free_dvector(ans,1,1);

  return bias;

}

void setHaloBiasIntegrand(double m, double y[], double deriv[], 
			double z1, Cosmology *c1, int flag)
{
  static double z;
  static Cosmology *c;
  double weight;

  if (flag == 1) {
    z = z1;
    c=c1;
    return;
  }

  weight=haloBias(m,z,c);
  deriv[1] = nm(m,z,c)*weight;
  
}

// Dummy function for integrating the bubble filling Q
void haloBiasIntegrand(double m, double y[], double deriv[])
{
  setHaloBiasIntegrand(m,y,deriv,0.0,NULL,0);
}

//number averaged bias for halos
double meanHaloBiasNW(double z, Cosmology *c, double mass)
{
    double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setHaloBiasIntegrandNW(0.0,ans,ans,z,c,1);
  mMin=coolMass(c,z);  //minimum halo mass

  if(mass>0.0) mMin=mass;

  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];

    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
    	   haloBiasIntegrandNW,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    mMin = mMax;
  }
  bias = ans[1];
  //bias /=c->nCollObject(z);  //normalise bias with number of collapsed objects
  bias/=meanHaloNumber(z,c);
  free_dvector(ans,1,1);

  return bias;

}

void setHaloBiasIntegrandNW(double m, double y[], double deriv[], 
			double z1, Cosmology *c1, int flag)
{
  static double z;
  static Cosmology *c;
  double weight;

  if (flag == 1) {
    z = z1;
    c=c1;
    return;
  }

  weight=haloBias(m,z,c);
  deriv[1] = nm(m,z,c)*weight/m;
  
}

// Dummy function for integrating the bubble filling Q
void haloBiasIntegrandNW(double m, double y[], double deriv[])
{
  setHaloBiasIntegrandNW(m,y,deriv,0.0,NULL,0);
}


//fit to bias from Tinker 2008
double biasTinker(double mass, double z, Cosmology *cosm)
{
   double nu, bias,dummy;
   double deltasc,sig;
   double A(1.0), a(0.1325);
   double B(0.183), b(1.5);
   double C(0.265), c(2.4);

   deltasc = cosm->delCrit0(z)/cosm->growthFac(z);
   sig = sigm(mass,dummy,cosm,0);
   nu=deltasc/sig;

   bias=1.0;
   bias-= A*pow(nu,a)/(pow(nu,a)+pow(deltasc,b));
   bias+=B*pow(nu,b);
   bias+=C*pow(nu,c);

   return bias;
}


//faster calculation of halo number
//mass averaged bias for halos
double meanHaloNumber(double z, Cosmology *c, double mass)
{
    double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setHaloNumberIntegrand(0.0,ans,ans,z,c,1);
  mMin=mass;
  if(mMin<0.0) mMin=coolMass(c,z);  //minimum halo mass


  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];

    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
    	   haloNumberIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    mMin = mMax;
  }
  bias = ans[1];
  //bias/=(CRITDENMSOLMPC*c->getOm0hh());
  //bias /= c->fColl(z);  //normalise halo mass function appropriately
  free_dvector(ans,1,1);

  return bias;

}

void setHaloNumberIntegrand(double m, double y[], double deriv[], 
			double z1, Cosmology *c1, int flag)
{
  static double z;
  static Cosmology *c;
  double weight;

  if (flag == 1) {
    z = z1;
    c=c1;
    return;
  }

  weight=1.0/m;
  deriv[1] = nm(m,z,c)*weight;
  
}

// Dummy function for integrating the bubble filling Q
void haloNumberIntegrand(double m, double y[], double deriv[])
{
  setHaloNumberIntegrand(m,y,deriv,0.0,NULL,0);
}
