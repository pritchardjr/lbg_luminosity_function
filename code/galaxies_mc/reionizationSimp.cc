//
// DISCLAIMER:  Currently my bubble distributions disagree with those
// of FZH2004 by ~20%.  The reason for this is still unclear.
//
//


/*  reionization.cc
 *
 *  contains code for reproducing the Furlanetto, Zaldarriaga, Hernquist
 *  bubble model for expanding HII regions.
 *
 *  contains other code useful for reionization
 *
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include "astrophysics.h"
#include "dcosmology.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "dnumrecipes.h"
#include "reionization.h"
#include "haloDensity.h"

using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

Reionization::Reionization(Cosmology *c1, Astrophysics *a1)
{
  double mb(1.0),dsdmB;
  // cout <<"Reionization constructor has been called." <<endl;
  c=c1;
  a=a1;
  sigm(mb,dsdmB,c,2);
  setZeta(a->getZeta());

  //  XICUT=1.0e-3;
  XICUT=1.0e-5;
  ASTROPHYSICS_FLAG=1;
 if(ASTROPHYSICS_FLAG==1) cout<<"normalizing to astrophysics ionization history"<<endl;
  //cout << "Reionization Constructor called successfully." <<endl; 
  
}

Reionization::~Reionization()
{
  // cout <<"Reionization destructor has been called." <<endl;

}

/////////////////////////////////////////////////////////////////////
// Member functions
/////////////////////////////////////////////////////////////////////
double Reionization::getZeta()
{
  return zetaBub;
}

// Set ionizing efficiency zetaBub
void Reionization::setZeta(double zetaIn)
{
  zetaBub=zetaIn;
}

double Reionization::getXICUT()
{
  return XICUT;
}

/////////////////////////////////////////////////////////////////////
// Fudge function
/////////////////////////////////////////////////////////////////////
//normalise zeta so that x_i from astrophysics class
// agrees with x_i from reionization model
void Reionization::normZeta(double z)
{
  double norm;
  double xiA, xiR;

  xiR=getZeta()*c->fColl(z);
  xiA=a->getXI(z);

  norm=xiA/xiR;

  zetaBub*=norm;
  cout<<getZeta()<<"\t"<<a->getZeta()<<"\t"<<norm<<endl;
}

/////////////////////////////////////////////////////////////////////
// Inverse error function code 
/////////////////////////////////////////////////////////////////////
//use GSL code to generate error function
//
//  given y=erf(x) find and return x=erf^{-1}(y)
//
double inverseERF(double y)
{
  double x;
  double x1(0.0),x2(10.0);
  double tol(1.0e-5);

  setIERF(y,y,1);

  zbracSimp(findIERF,x1,x2);
  x=zriddrSimp(findIERF,x1,x2,tol);
  
  return x;
}

double findIERF(double x)
{
  return setIERF(x,0.0,0);
}

double setIERF(double x, double y1, int iflag)
{
  static double y;

  if(iflag==1){
    y=y1;
    return 0.0;
  }

  return gsl_sf_erf(x)-y;
}

/////////////////////////////////////////////////////////////////////
// Bubble distribution code - FZH2004
/////////////////////////////////////////////////////////////////////
// Calculate the mass distribution of bubbles
// dndlm - number density per log interval of mass
//
// Units: Mpc^{-3} (Comoving)
double Reionization::nmBub(double mb, double z)
{
  double K,deltac;
  double B,B0,B1;
  double sigMbub,sigMmin;
  double dsdmB,dsdmMin;
  double dndlm;
  double rho;
  double mMin;
  double Q;
  static double renorm(1.0);
  static int count;
  static double zsave(-1.0);
  static double zetaSave(-1.0);

  //renormalise output to ensure Q=zeta*fcoll
  if(count==0 || fabs(z-zsave)>1.0e-4 || fabs(zetaBub-zetaSave)>1.0e-4){
    zsave=z;
    if(ASTROPHYSICS_FLAG==1){
      normZeta(z);  // normalise zeta to match ionization history
    }
    zetaSave=zetaBub;
    count++;
    renorm=1.0;
    Q=fillingQ(z);
    renorm=getZeta()*c->fColl(z)/Q;
    cout<<"renorm \t"<<z<<"\t"<<Q<<"\t"<<getZeta()*c->fColl(z)<<"\t"<<renorm<<endl;
  }


  deltac=c->delCrit0(z)/c->growthFac(z);
  K=inverseERF(1.0-1.0/getZeta());
  mMin=coolMass(c,z);  //minimum mass of luminous object

  sigMbub=sigm(mb,dsdmB,c,1);
  sigMmin=sigm(mMin,dsdmMin,c,1);
  
  B0=deltac-sqrt(2.0)*K*sigMmin;
  if(B0<1.0e-20){
    // cout<<"Fully ionized"<<endl;
    return 0.0;  //at low z B0 can turn negative indicating reion
  }
  B1=K/(sqrt(2.0)*sigMmin);
  B=B0+B1*sigMbub*sigMbub;

  rho=CRITDENMSOLMPC*c->getOm0hh();

  dndlm=sqrt(2.0/PI)*rho/mb;
  dndlm*=fabs(dsdmB*mb/sigMbub);
  dndlm*=B0/sigMbub;
  dndlm*=exp(-B*B/2.0/sigMbub/sigMbub);

  // cout<<z<<"\t"<<mb<<"\t"<<B0<<"\t"<<B1<<"\t"<<sigMbub<<"\t"<<dndlm<<endl;

  dndlm*=renorm;

  return dndlm;
}

//mean filling fraction of bubbles \bar{Q}
//
double Reionization::fillingQ(double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double Q;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setQIntegrand(0.0,ans,ans,z,this,1);
  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble
  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
	   qIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6 ||fabs(ans[1])<1.0e-60)
      break;
    mMin = mMax;
  }
  Q = ans[1]/(CRITDENMSOLMPC*c->getOm0hh());
  free_dvector(ans,1,1);

  if(Q>1.0||Q<0.0){
    cout<<"Bubbles fill all space"<<endl;
    Q=1.0;
  }

  return Q;
}

// Integrand for computing the bubble filling fraction.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setQIntegrand(double m, double y[], double deriv[], 
			double z1, Reionization *rz1, int flag)
{
  static double z;
  static Reionization *rz;
  double qDensity;

  if (flag == 1) {
    z = z1;
    rz = rz1;
    return;
  }

  //qDensity = rz->nmBub(m,z)/m;
  //deriv[1] = qDensity*m;

  deriv[1]=rz->nmBub(m,z);

}

// Dummy function for integrating the bubble filling Q
void qIntegrand(double m, double y[], double deriv[])
{
  setQIntegrand(m,y,deriv,0.0,NULL,0);
}

//mean bubble size <R>
//
double Reionization::meanRBub(double z)
{
  double Rbub;
  double muse,mMin,mMax;
  double tol(1.0e-4);
  double step;
  int i, nstep(3);
  double nBub;

  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble
  mMax=1.0e10*mMin;

  setRBub(1.0,z,this,c,1);
  step=exp(log(0.99999*mMax/mMin)/(double)(nstep));
  Rbub=0.0;
  muse=mMin;
  for(i=1;i<=nstep;i++){
    Rbub+=qromb(intRBub,log(muse),log(muse*step),tol); 
    muse*=step;
  }
  
  nBub=meanNBub(z);
  //  cout<<Rbub<<"\t"<<nBub<<endl;
  Rbub/=nBub;
  return Rbub;
}

double setRBub(double lmb, double z1, Reionization *rz1, Cosmology *c1, int iflag)
{
  static double z;
  static Reionization *rz;
  static Cosmology *c;
  double dndlm;
  double rho;
  double dRdlm;
  double mb;
  double R,V;

  if(iflag==1){
    z=z1;
    c=c1;
    rz=rz1;
    return 0.0;
  }

  mb=exp(lmb);

  dndlm=rz->nmBub(mb,z);
  rho=CRITDENMSOLMPC*c->getOm0hh();
  V=mb/rho;
  R=pow(V/(4.0*PI/3.0),0.33333);
  
  dRdlm=dndlm*R;  //This needs to be normalised correctly by 
                  // 1/ (int dndlm)
  //  cout<<mb<<"\t"<<R<<"\t"<<dndlm<<endl;

  if(iflag==2){
    dRdlm=dndlm;
  }
  
  return dRdlm*V;
}

double intRBub(double lmb)
{
  return setRBub(lmb,0.0,NULL,NULL,0);
}

//integrate the bubble number density
double Reionization::meanNBub(double z)
{
  double Nbub;
  double muse,mMin,mMax;
  double tol(1.0e-4);
  double step;
  int i, nstep(3);

  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble
  mMax=1.0e10*mMin;

  setRBub(1.0,z,this,c,1);
  step=exp(log(0.99999*mMax/mMin)/(double)(nstep));
  Nbub=0.0;
  muse=mMin;
  for(i=1;i<=nstep;i++){
    Nbub+=qromb(intNBub,log(muse),log(muse*step),tol); 
    muse*=step;
  }
  return Nbub;
}

double intNBub(double lmb)
{
  return setRBub(lmb,0.0,NULL,NULL,2);
}

//calculate barrier for bubble model
double Reionization::barrierBub(double mb, double z, double BB[])
{
  double K,deltac;
  double B,B0,B1;
  double sigMbub,sigMmin;
  double dsdmB,dsdmMin;
  double mMin;

  deltac=c->delCrit0(z)/c->growthFac(z);
  K=inverseERF(1.0-1.0/getZeta());
  mMin=coolMass(c,z);  //minimum mass of luminous object

  sigMbub=sigm(mb,dsdmB,c,1);
  sigMmin=sigm(mMin,dsdmMin,c,1);
  
  B0=deltac-sqrt(2.0)*K*sigMmin;
  if(B0<1.0e-5){
    cout<<"Fully ionized"<<endl;
    return 0.0;  //at low z B0 can turn negative indicating reion
  }
  B1=K/(sqrt(2.0)*sigMmin);
  B=B0+B1*sigMbub*sigMbub;

  BB[1]=B0;
  BB[2]=B1;
  BB[3]=B;

  return B*c->growthFac(z);  //physical overdensity
}


// Calculate the bubble bias for bubbles of mass mb at redshift z
//
//
double Reionization::biasBub(double mb, double z)
{
  double bias;
  double *BB;
  double B0,B;
  double sigMbub,dsdmB;
  BB=dvector(1,3);

  barrierBub(mb,z,BB);
  B0=BB[1];
  B=BB[3];
  free_dvector(BB,1,3);

  sigMbub=sigm(mb,dsdmB,c,1);

  bias=1.0+(B/sigMbub/sigMbub-1.0/B0)/c->growthFac(z);
  // cout<<mb<<"\t"<<z<<"\t"<<bias<<"\t"<<BB[2]*sigMbub*sigMbub<<"\t"<<B0<<endl;
  // bias=1.0;
  return bias;
}

// Calculate the bubble bias for bubbles of mass mb at redshift z
//
//  THIS CURRENTLY USES THE UNCORRECTED EXPRESSION FROM FZH04
// keep it around for comparison purposes
double Reionization::biasBubFZH(double mb, double z)
{
  double bias;
  double *BB;
  double B0,B;
  double sigMbub,dsdmB;
  BB=dvector(1,3);

  barrierBub(mb,z,BB);
  B0=BB[1];
  B=BB[3];
  free_dvector(BB,1,3);

  sigMbub=sigm(mb,dsdmB,c,1);

  bias=1.0+B0*B0/sigMbub/sigMbub/B;

  return bias;
}


//mean bias of bubbles
//
double Reionization::meanBiasBub(double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double Q;
  double bias;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setBiasBubIntegrand(0.0,ans,ans,z,this,1);
  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble
  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
	   biasBubIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6)
      break;
    mMin = mMax;
  }
  bias = ans[1]/(CRITDENMSOLMPC*c->getOm0hh());
  Q=fillingQ(z);
  bias/=Q;
  free_dvector(ans,1,1);

  return bias;
}

// Integrand for computing the bubble filling fraction.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setBiasBubIntegrand(double m, double y[], double deriv[], 
			double z1, Reionization *rz1, int flag)
{
  static double z;
  static Reionization *rz;
  double qDensity;

  if (flag == 1) {
    z = z1;
    rz = rz1;
    return;
  }

  // qDensity = rz->nmBub(m,z)/m;
  //deriv[1] = qDensity*m*rz->biasBub(m,z);
  deriv[1] = rz->nmBub(m,z)*rz->biasBub(m,z);

}

// Dummy function for integrating the bubble filling Q
void biasBubIntegrand(double m, double y[], double deriv[])
{
  setBiasBubIntegrand(m,y,deriv,0.0,NULL,0);
}



/////////////////////////////////////////////////////////////////////
// FZH04 HII region correlation functions
/////////////////////////////////////////////////////////////////////
//Calculate full Brightness temperature correlation function
double Reionization::correlateTb(double r, double z)
{
  double corrXX, corrDD, corrXD, corr;
  double xH;

  xH=1.0-fillingQ(z);
 
  corrXX=correlateXX(r,z);
  corrXD=correlateXD(r,z);
  corrDD=xidd(r,z,c,0);
 
  corr=corrXX*(1.0+corrDD);
  corr+= xH*xH*corrDD;
  corr+= corrXD*(2.0*xH+corrXD);
  
  return corr;
}




//Calculate the xx correlation function

//Calculate overlap volume
double Reionization::overlapVolume(double r, double mb)
{
  double vb,rb;
  double vo;

  vb=mb/(CRITDENMSOLMPC*c->getOm0hh());
  rb=pow(vb/(4.0*PI/3.0),1.0/3.0);

  if(r>2.0*rb){
    vo=0.0;
  }else{
    vo=vb-PI*r*(rb*rb-r*r/12.0);
  }
  //cout<<rb<<"\t"<<r<<"\t"<<vb<<"\t"<<vo<<endl;
  return vo;
}

//Calculate full XX correlation function
double Reionization::correlateXX(double r, double z)
{
  double corr1b, corr2b,corrXX,Q;
  double renorm;
  double overlap,overlapb;
  double overlapd;
  double exparg,temp;
  double biasb;
  double highXX;
  static double rbreak(1.0e30);
  static double zsave(-1.0);
  
  //renormalise to correct for overlap
  Q=fillingQ(z);
  if(Q<XICUT) return 0.0;  //no bubbles
  if(fabs(Q-1.0)<1.0e-2) return 0.0;  //fully ionized
  if(fabs(z-zsave>1.0e-4)){
    rbreak=1.0e30;
    zsave=z;
  }  

  renorm= -log(1.0-Q)/Q;
  biasb=meanBiasBub(z);
  highXX=(1.0-Q)*renorm*pow(Q*biasb,2.0)*xidd(r,z,c,0);

  if(r>rbreak) return highXX;

  //  cout<<"overlap...";
  overlap=meanBubOverlap(r,z);
  //  cout<<"overlap_bias...";

  overlapb=biasb*xidd(r,z,c,0)*(Q*biasb-meanBubOverlapBias(r,z));
  overlapd=Q-overlap;

  //c loses numerical precision on 1-exp(-x) once x<1,0e-16
  //so need to be a little careful on how we handle these quantities
  //for x<1.0e-7 can use approx 1-exp(-x)=x without obvious loss of precision


  exparg=renorm*overlap;
  if(fabs(exparg)<1.0e-8){
    corr1b=exparg;
  }else{
    corr1b=1.0-exp(-exparg);
  }

  corr2b=exp(-renorm*overlap);
  exparg=renorm*overlapd;
  if(fabs(exparg)<1.0e-8){
    corr2b*=exparg;
  }else{
    corr2b*=(1.0-exp(-exparg));
  }

  if(fabs(renorm*overlapb)<1.0e-8 && fabs(renorm*(Q-overlapd))<1.0e-8){
    corr2b*=Q+(1.0-Q)*renorm*overlapb;
  }else{
    corr2b*=1.0-exp(-renorm*overlapd)*exp(-renorm*overlapb);
  }

  corrXX=corr1b+corr2b-Q*Q;
  //check to see if we can apply limiting behaviour 
  if(fabs((corrXX-highXX)/(corrXX))<1.0e-6){
    rbreak=r;
  }
  
  //  cout<<"overlap:"<<overlap<<"\t"<<"Q- overlap:"<<overlapd<<"\t"<<"bias overlap:"<<overlapb<<endl;
  //  cout<<"XX:"<<corr1b+corr2b-Q*Q<<"\t"<<Q<<"\t"<<corr1b<<"\t"<<corr2b<<"\t"<<corr2b/Q<<"\t"<<(corr1b+corr2b-Q*Q)/highXX<<endl;

  return corrXX;
}

//Calculate XX correlation function where two points separated by 
//a distance r are ionized by the same source (1 bubble term)
double Reionization::correlateXX1b(double r, double z)
{
  double corr;
  double renorm,Q;
  double overlap;

  //renormalise to correct for overlap
  Q=fillingQ(z);
  renorm= -log(1.0-Q)/Q;

  overlap=meanBubOverlap(r,z);

  corr=1.0-exp(-renorm*overlap);
  return corr;
}

//Calculate XX correlation function where two points separated by 
//a distance r are ionized by the different sources (2 bubble term)
double Reionization::correlateXX2b(double r, double z)
{
  double corr;
  double renorm,Q;
  double exparg;

  //renormalise to correct for overlap
  Q=fillingQ(z);
  renorm= -log(1.0-Q)/Q;

  corr=exp(-renorm*meanBubOverlap(r,z));
  exparg=renorm*(fillingQ(z)-meanBubOverlap(r,z));

  //c loses numerical precision on 1-exp(-x) once x<1,0e-16
  //so need to be a little careful on how we handle these quantities
  //for x<1.0e-7 can use approx 1-exp(-x)=x without obvious loss of precision
  if(fabs(exparg)<1.0e-8){
    corr*=exparg;
  }else{
    corr*=(1.0-exp(-exparg));
  }

  exparg=renorm*meanBubOverlapCluster(r,z);
  if(fabs(exparg)<1.0e-8){
    corr*=exparg;
  }else{
    corr*=(1.0-exp(-exparg));
  }

  return corr;
}

//Calculate X-density cross-correlation function
double Reionization::correlateXD(double r, double z)
{
  double corr(0.0);
  double Q;
  double prob;
  double highXD,lowXD;
  double halobias;
  double meanbias;
  static double rbreak(1.0e30);
  static double zsave(-1.0);
  
  Q=fillingQ(z);
  if(Q<XICUT) return 0.0;  //no bubbles
  if(fabs(Q-1.0)<1.0e-2) return 0.0;  //fully ionized
  if(fabs(z-zsave>1.0e-4)){
    rbreak=1.0e30;
    zsave=z;
  }

  halobias=meanHaloBias(z,c);
  meanbias=meanBiasBub(z);
  highXD=(1.0-Q)*meanbias*Q*xidd(r,z,c,0)*halobias;

  if(r>rbreak) return -highXD/halobias;
  //  lowXD=(1.0-Q)*(1.0-exp(-meanbias*Q*xidd(meanRBub(z),z,c,0)*halobias));
  
  //  cout<<"integrating...";
  prob=correlateXDIntegral(r,z);
  corr=(1.0-Q)*prob;

  //   cout<<"XD:"<<r<<"\t"<<-corr/halobias<<"\t"<<corr<<"\t"<<lowXD<<"\t"<<highXD<<"\t"<<corr/highXD<<"\t"<<corr/lowXD<<endl;

  //once reach numerical average can use approximation to spead things up
  if(fabs((corr-highXD)/corr)<1.5e-8){
    rbreak=r;
  }

  return -corr/halobias;  //we want neutral-density not ionized-density
}


//mean overlap of bubbles
//
double Reionization::meanBubOverlap(double r, double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;
  double rMin;
  double rb;
  //  static int count;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setBubOverlapIntegrand(0.0,ans,ans,z,r,this,1);
  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble

  //integrator will complain if kernel is zero for long periods
  //If we're going to fall in a regime where the smallest bubbles are still
  //too large to get overlap increase the start of integration to the smallest
  //bubble size where overlap occurs, then increase by small fudge fraction
  rMin=pow(mMin/(CRITDENMSOLMPC*c->getOm0hh())/(4.0*PI/3.0),1.0/3.0);
  if(r>2.0*rMin){
    mMin=4.0*PI/3.0*CRITDENMSOLMPC*c->getOm0hh()*pow((r+1.0e-4)/2.0,3.0);
    //    cout<<"changing mMin in overlap"<<endl;
  }
   
  //check to see if we've reached regime where nm_bub=0
  // if(nmBub(mMin,z)<1.0e-20) return 0.0; 

  //count=0;
  while (1) {
    //count++;
    //cout<<count<<".";
    step = mMin/10.0;
    mMax = mMin*10.0;

    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
	   bubOverlapIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6 || fabs(ans[1])<1.0e-50){
      //cout<<"break";
      break;
    }   
    mMin = mMax;

  }
  //  cout<<endl;
  bias = ans[1];
  free_dvector(ans,1,1);

  return bias;
}

// Integrand for computing the bubble filling fraction.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setBubOverlapIntegrand(double m, double y[], double deriv[], 
			double z1, double r1, Reionization *rz1, int flag)
{
  static double z;
  static double r;
  static Reionization *rz;

  if (flag == 1) {
    z = z1;
    rz = rz1;
    r=r1;
    return;
  }

  deriv[1] = rz->nmBub(m,z)/m*rz->overlapVolume(r,m);
  //  cout<<"\t"<<m<<"\t"<<rz->overlapVolume(r,m)<<"\t"<<deriv[1]<<"\t"<< rz->nmBub(m,z)<<endl;

}

// Dummy function for integrating the bubble filling Q
void bubOverlapIntegrand(double m, double y[], double deriv[])
{
  setBubOverlapIntegrand(m,y,deriv,0.0,0.0,NULL,0);
}

//mean overlap of bubbles with Clustering
//
double Reionization::meanBubOverlapCluster(double r, double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;
  //  int count;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setBubOverlapClusterIntegrand(0.0,ans,ans,z,r,c,this,1);
  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble

  //count=0;
  while (1) {
    //count++;
    //cout<<count<<".";
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    //  if(nmBub(mMax,z)<1.0e-40) break;
    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
	   bubOverlapClusterIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4 || fabs(ans[1])<1.0e-50) break;
    mMin = mMax;
  }

  
  bias = ans[1];
  free_dvector(ans,1,1);
  //  cout<<endl;

  return bias;
}

// Integrand for computing the bubble filling fraction.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setBubOverlapClusterIntegrand(double m, double y[], double deriv[], 
			double z1, double r1, Cosmology *c1,Reionization *rz1, int flag)
{
  static double z;
  static double r;
  static Reionization *rz;
  static Cosmology *c;
  double vb,xi,vo;

  if (flag == 1) {
    z = z1;
    rz = rz1;
    c=c1;
    r=r1;
    return;
  }
  vb=m/(CRITDENMSOLMPC*c->getOm0hh());
  vo=rz->overlapVolume(r,m);
  xi=rz->correlateMeanBub(r,m,z);
  deriv[1] = rz->nmBub(m,z)/m*(vb-vo)*(1.0+xi);
  // cout<<m<<"\t"<<rz->nmBub(m,z)<<"\t"<<xi<<"\t"<<vb<<"\t"<<vo<<"\t"<<deriv[1]<<endl;

}

// Dummy function for integrating the bubble filling Q
void bubOverlapClusterIntegrand(double m, double y[], double deriv[])
{
  setBubOverlapClusterIntegrand(m,y,deriv,0.0,0.0,NULL,NULL,0);
}

//Calculate the excess probability that the point r2 is ionized
//by a source of size m, given that r1 is ionized.  r=r2-r1
//  FZH04:  \bar{\zeta}
//

double Reionization::correlateMeanBub(double r, double mb, double z)
{
  static int count;
  static double zsave(-1.0);
  double correlate;
  double bias, xi;
  static double renorm;
  double Q;

  //initialise xidd lookup table if change z or when starting
  if (count==0 || fabs(z-zsave)>1.0e-4){
    zsave=z;
    count++;
    cout<<"Initialising xidd lookup table"<<endl;
    xidd(r,z,c,1); 
    cout<<"xidd lookup initialised"<<endl;
    Q=fillingQ(z);
    renorm= -log(1.0-Q)/Q;
  }

  xi=xidd(r,z,c,0);   //use haloDensity code for this
  bias=meanBiasBub(z)*renorm;  //factor of nmBub in meanBiasBub
  correlate=bias*biasBub(mb,z)*xi;

  return correlate;
}



//Calculate double integral for XD cross-correlation
//
// Integral over halo mass not bubble mass
//
double Reionization::correlateXDIntegral(double r, double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;
  //static int count;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setCorrXDIntegrand(0.0,ans,ans,z,r,c,this,1);
  mMin=coolMass(c,z);  //minimum halo mass

  //  count=0;
  while (1) {
    //    count++;
    //cout<<count<<".";
    step = mMin/10.0;
    mMax = mMin*10.0;

    //  cout<<mMin<<"\t"<<mMax<<"\t"<<oldans<<"\t"<<ans[1]<<endl;

    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
    	   corrXDIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6  || fabs(ans[1])<1.0e-50){
      //cout<<"break";
      break;
    }
    mMin = mMax;
  }
  bias = ans[1]/(CRITDENMSOLMPC*c->getOm0hh());
  bias /= c->fColl(z);  //normalise halo mass function appropriately
  free_dvector(ans,1,1);
  //cout<<endl;

  return bias;
}

// Integrand for computing the bubble filling fraction.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setCorrXDIntegrand(double m, double y[], double deriv[], 
			double z1, double r1, Cosmology *c1,Reionization *rz1, int flag)
{
  static double z;
  static double r;
  static Reionization *rz;
  static Cosmology *c;
  static double Q;
  double weight;
  static double renorm;
  static double xiAverage;
  static double bubbias;
  double exparg;

  double epa;

  if (flag == 1) {
    z = z1;
    c=c1;
    rz = rz1;
    r=r1;
    Q=rz->fillingQ(z);
    //renormalise to correct for overlap 
    renorm= -log(1.0-Q)/Q;
    renorm=1.0;
    // xiAverage is mass independent
    //note also that xiAverage contains all r dependence
    //    cout<<"XiXD...";
    xiAverage=rz->meanXiXD(r,z);
    //epa=Q*rz->meanBiasBub(z)*xidd(r,z,c,0);
    //cout<<endl;
    //cout<<"XDcheck \t"<<epa<<"\t"<<xiAverage<<"\t"<<xiAverage/epa<<endl;
    //cout<<"got..."<<xiAverage;
    return;
  }

  exparg=renorm*haloBias(m,z,c)*xiAverage;

  if(fabs(exparg)<1.0e-8){
    weight=exparg;
  }else {
    weight=1.0-exp(-exparg);
  }
  //  cout<<"w:"<<weight<<"\t"<<renorm<<endl;
  deriv[1] = nm(m,z,c)*weight;
  
}

// Dummy function for integrating the bubble filling Q
void corrXDIntegrand(double m, double y[], double deriv[])
{
  setCorrXDIntegrand(m,y,deriv,0.0,0.0,NULL,NULL,0);
}

//Calculate <\xi> Equation (24) from FZH04
//
double Reionization::meanXiXD(double r, double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setXiXDIntegrand(0.0,ans,ans,z,r,c,this,1);
  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble
  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    //  cout<<mMin<<"\t"<<mMax<<"\t"<<oldans<<"\t"<<ans[1]<<"\t"<<fabs((ans[1]-oldans)/ans[1])<<endl;
    oldans = ans[1];
 
    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
	   xiXDIntegrand,bsstep);
    if ((fabs((ans[1]-oldans)/ans[1]) < 1.0e-6  || fabs(ans[1])<1.0e-50))
      break;
    // cout<<ans[1]<<"\t"<<oldans<<"\t"<<(ans[1]-oldans)/ans[1]<<endl;
    mMin = mMax;
  }
  //  cout<<mMin<<"\t"<<mMax<<"\t"<<oldans<<"\t"<<ans[1]<<"\t"<<fabs((ans[1]-oldans)/ans[1])<<endl;
  bias = ans[1];
  free_dvector(ans,1,1);

  return bias;
}

// Integrand for computing <\xi>.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setXiXDIntegrand(double m, double y[], double deriv[], 
			double z1, double r1, Cosmology *c1, Reionization *rz1, int flag)
{
  static double z;
  static double r;
  static Reionization *rz;
  static Cosmology *c;
  double vb,rb,reff,xi;

  if (flag == 1) {
    z = z1;
    rz = rz1;
    c=c1;
    r=r1;
    return;
  }

  vb=m/(CRITDENMSOLMPC*c->getOm0hh());
  rb=pow(vb/(4.0*PI/3.0),1.0/3.0); 
  reff=fmax(rb,r);
  xi=xidd(reff,z,c,0)*vb;

  //   deriv[1] =  rz->nmBub(m,z)/m*rz->biasBub(m,z)*rz->xiAverageInt(m,r,z);
        deriv[1] =  rz->nmBub(m,z)/m*rz->biasBub(m,z)*xi; 
   //  cout<<m<<"\t"<<rz->nmBub(m,z)/m<<"\t"<<deriv[1]<<"\t"<<rz->xiAverageInt(m,r,z)<<"\t"<<rz->biasBub(m,z)<<endl;

}

// Dummy function for integrating the <\xi>
void xiXDIntegrand(double m, double y[], double deriv[])
{
  setXiXDIntegrand(m,y,deriv,0.0,0.0,NULL,NULL,0);
}

 double rsaveXD;
double rXD;
//int xdcount;
//int xdcountr;

// Calculate integral over DD correlation function over volume of a bubble
// of mass mb
//
double Reionization::xiAverageInt(double mb, double r, double z)
{
  double vb,rb;
  double rMin,rMax;
  double tol(1.0e-4);
  double result;

  //limits of integration
  vb=mb/(CRITDENMSOLMPC*c->getOm0hh());
  rb=pow(vb/(4.0*PI/3.0),1.0/3.0); 

  rXD=r;

  //   xdcount=0;
  //   xdcountr=0;
  setXiAverageIntegrand(0.0,0.0,r,z,c,1);
  rMax=rb;
  rMin=1.0e-5;  
  //  cout<<"here"<<r<<"\t"<<rMin<<"\t"<<rMax<<"\t"<<(r-rMax)<<"\t"<<(r+rMax)<<endl;
  result=qromb(intXiAverageR,rMin,rMax,tol);

  //  cout<<"there"<<"\t"<<xdcount<<"\t"<<xdcountr<<"\t"<<xdcount/xdcountr<<endl;
  result*=2.0*PI;

  return result;
}


double intXiAverageR(double rs)
{
  double tol(1.0e-4);
  double temp;
  double cutoff(1.0),mucut(0.95);
   int waymark,start;
  rsaveXD=rs;
  //   xdcountr++;
  
  if(fabs(rXD-rs)<cutoff){
    //  cout<<"splitting integral"<<endl;
    // start=xdcount;
    temp= qromb(intXiAverageMu,-1.0,mucut,tol);
    //  waymark=xdcount;
    temp+= qromb(intXiAverageMu,mucut,1.0,tol);
    //  cout<<"\t"<<waymark<<"\t"<<waymark-start<<"\t"<<xdcount-waymark<<endl;
  }else{
    temp= qromb(intXiAverageMu,-1.0,1.0,tol);
  }
 
   return temp;
}

double intXiAverageMu(double mu)
{
  return setXiAverageIntegrand(rsaveXD,mu,0.0,0.0,NULL,0);
}

double setXiAverageIntegrand(double rs, double mu, double r1, double z1, Cosmology *c1,int iflag)
{
  static double r;
  static double z;
  static Cosmology *c;
  double d,xi;

  if(iflag==1)
    {
      r=r1;
      z=z1;
      c=c1;
      return 0.0;
    }

  d=sqrt(rs*rs+r*r-2.0*rs*r*mu);
  xi=xidd(d,z,c,0);

  return xi*rs*rs;
}

//////////////////////////////////
double musaveXD, rminXD, rmaxXD;
double rpivotXD;

double Reionization::xiAverageInt2(double mb, double r, double z)
{
  double vb,rb;
  double rMin,rMax;
  double tol(1.0e-4);
  double result;

  //limits of integration
  vb=mb/(CRITDENMSOLMPC*c->getOm0hh());
  rb=pow(vb/(4.0*PI/3.0),1.0/3.0); 

  setXiAverageIntegrand(0.0,0.0,r,z,c,1);
  rmaxXD=rb;
  rminXD=1.0e-5;  
  rpivotXD=r;
  
  //  if(rmaxXD>r)cout<<"difficult"<<endl;
  result=qromb(intXiAverageMu2,-1.0,1.0,tol);

  result*=2.0*PI;

  return result;
}


double intXiAverageR2(double rs)
{
  return setXiAverageIntegrand(rs,musaveXD,0.0,0.0,NULL,0);
}

double intXiAverageMu2(double mu)
{
  double temp, tol(1.0e-4);

  musaveXD=mu;
  temp=qromb(intXiAverageR2,rminXD,rmaxXD,tol); 
  //  temp=qromb(intXiAverageR2,rpivotXD*mu,rminXD,tol);
  //temp+=qromb(intXiAverageR2,rpivotXD*mu,rmaxXD,tol);

  return temp;
}


///////////////////////////////////////////////////////////////////
//  McQuinn HII correlation model
///////////////////////////////////////////////////////////////////
//
// McQuinn et al 2005 produced an improved model for the HII 
// correlation functions and corrected a few small mistakes in FZH04.
// This code reproduces their results

//NOT FINISHED
double Reionization::correlateXXMQ(double r, double z)
{
  double P1, P2, Q;
  double renorm;
  double eta(0.5); //transition neutral fraction for approximations
  double xi;

  Q=fillingQ(z);
  renorm=-log(1.0-Q)/Q;
  P1=renorm*meanBubOverlap(r,z);

  if(Q>eta){
    xi=(1.0-Q)*P1+Q*Q;
  }else{
    P2=renorm*renorm*p2MQ(r,z);
    xi=P1+P2;
  }
  return xi;
}

double Reionization::correlateXDMQ(double r, double z)
{
  double P1, Pout,Pin, Q;
  double renorm, overlap;
  double eta(0.5); //transition neutral fraction for approximations
  double xi;

  Q=fillingQ(z);
  renorm=-log(1.0-Q)/Q;
  
  overlap=meanBubOverlap(r,z);
  Pin=renorm*(overlap+meanBubOverlapBias(r,z));

  if(Q>eta){
    P1=renorm*overlap;
    xi=Pin-P1;
  }else{
    Pout=Q-renorm*overlap+renorm*meanXiXD(r,z);
    xi=Pin+Pout-Q;
  }
  return xi;
}


double Reionization::p2MQ(double r, double z)
{
  return 1.0;
}

//mean overlap of bubbles with extra bias factor
//
double Reionization::meanBubOverlapBias(double r, double z)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double bias;
  double rMin;

  ans = dvector(1,1);
  ans[1] = 0.0;
  setBubOverlapBiasIntegrand(0.0,ans,ans,z,r,c,this,1);
  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble

  //integrator will complain if kernel is zero for long periods
  //If we're going to fall in a regime where the smallest bubbles are still
  //too large to get overlap increase the start of integration to the smallest
  //bubble size where overlap occurs, then increase by small fudge fraction
  // 1.01 here.
  rMin=pow(mMin/(CRITDENMSOLMPC*c->getOm0hh())/(4.0*PI/3.0),1.0/3.0);
  if(r>=2.0*rMin){ 
    // mMin=1.01*4.0*PI/3.0*CRITDENMSOLMPC*c->getOm0hh()*pow(r/2.0,3.0);
    mMin=4.0*PI/3.0*CRITDENMSOLMPC*c->getOm0hh()*pow((r+1.0e-4)/2.0,3.0);
    //    cout<<"...changing mMin in biasOverlap"<<endl;
  }

  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-6,step,0.0,&goodSteps,&badSteps,
	   bubOverlapBiasIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-6 || fabs(ans[1])<1.0e-60){
      // cout<<"\t checking: "<<ans[1];
      break;
    }
    mMin = mMax;
    //  cout<<"\t checking: "<<ans[1]<<endl;;
  }
  bias = ans[1];
  free_dvector(ans,1,1);

  return bias;
}

// Integrand for computing the bubble filling fraction.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setBubOverlapBiasIntegrand(double m, double y[], double deriv[], 
			double z1, double r1, Cosmology *c1,Reionization *rz1, int flag)
{
  static double z;
  static double r;
  static Reionization *rz;
  static Cosmology *c;
  double vo;
  static double bias;

  if (flag == 1) {
    z = z1;
    rz = rz1;
    r=r1;
    c=c1;
    bias=rz->meanBiasBub(z);
    return;
  }
  vo=rz->overlapVolume(r,m);
 
  //using bubble bias causes problems at low x_H, where many decades
  //of bubble mass violate the assumption under which the MQ bias was
  //derived.  Replacing with the mean bubble bias gives a very similar
  // xx power spectrum (with increasing difference at low x_H).
  //  deriv[1] = rz->nmBub(m,z)/m*vo*rz->biasBub(m,z);
  deriv[1] = rz->nmBub(m,z)/m*vo*bias;

}

// Dummy function for integrating the bubble filling Q
void bubOverlapBiasIntegrand(double m, double y[], double deriv[])
{
  setBubOverlapBiasIntegrand(m,y,deriv,0.0,0.0,NULL,NULL,0);
}


/////////////////////////////////////////////////////////////////
// Improved 3D Fourier transform code
/////////////////////////////////////////////////////////////////

const double INT_ACCURACY2(1.0e-5);

int nrJ(100),nkJ(100);
double *rftJ,*xiftJ,*xiftspJ;

// Function to construct a table of the (3D, isotropic) FT of the
// function input as xivec.  
//   xivec = correlation function values
//   rvec = ordinates of the correlation function (comoving Mpc).
//   nrJ = number of points in xivec/rvec
//   pvec = vector for storing power spectrum
//   kvec = vector containing wavenumbers at which to evaluate P(k)
//          (comoving Mpc^-1)
//   nkJ = number of points in kvec/pvec
// Uses spline interpolation on the correlation function to do the 
// evaluation.  Integrates over setFTQRInt.
//
void FTTableJP(double xivec[], double rvec[], int nrJ1, double pvec[], 
	     double kvec[], int nkJ1)
{
  double rmin,rmax;
  double k(1.0);
  int ift;
  const double natural(1.0e30);
  double tol(1.0e-4);
  double L;

  nrJ=nrJ1;
  nkJ=nkJ1;

  // Set up the spline fit to the correlation function
  rftJ = dvector(1,nrJ);
  xiftJ = dvector(1,nrJ);
  xiftspJ = dvector(1,nrJ);

  for (ift=1;ift<=nrJ;ift++) {
    rftJ[ift] = rvec[ift];
    xiftJ[ift] = xivec[ift]*rvec[ift];  //include factor of r in spline
    // xiftJ[ift] = xivec[ift];
  }

  spline(rftJ,xiftJ,nrJ,natural,natural,xiftspJ);
  cout<<"spline OK"<<endl;

  rmin=rvec[1];
  rmax=rvec[nrJ];

  /////////////////
  double r, xi;
  ofstream fout;
  char *file;

  file="corrcheck.dat";
  fout.open(file);
  r=rmin;
  while(r<=rmax){
    splint(rftJ,xiftJ,xiftspJ,nrJ,r,&xi);
    //  cout<<r<<"\t"<<xi<<endl;
    fout<<r<<"\t"<<xi<<endl;
    r*=1.1;
  }
  fout.close();
  //////////////////

  
  rmin=rvec[1];
  rmax=rvec[nrJ];
  L=rmax-rmin;


  ///////////////////////////////////////////

  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc(1000);

  gsl_integration_workspace * cw 
    = gsl_integration_workspace_alloc(1000);

 gsl_integration_qawo_table * wf
   = gsl_integration_qawo_table_alloc (k,1.0,GSL_INTEG_SINE,1000); 
  
  double result, error;
  double alpha(1.0);

  gsl_function F;
  F.function = &FTQRIntJP;
  F.params = &alpha;


  for (ift=1;ift<=nkJ;ift++) {
    cout<<ift<<".";
    k=kvec[ift];
    setFTQRIntJP(0.0,k,rvec[1],rvec[nrJ],1);
    //integration routine
    gsl_integration_qawo_table_set (wf,k,L,GSL_INTEG_SINE);
    gsl_integration_qawo (&F,rmin,0,tol, 1000,
                        w, wf,&result, &error); 
    pvec[ift] = 4.0*PI*result;
  }
  cout<<"... FT error="<<error<<endl;

  free_dvector(rftJ,1,nrJ);
  free_dvector(xiftJ,1,nrJ);
  free_dvector(xiftspJ,1,nrJ);
  return;
}


// Integrand for making the power spectrum, just a 3D isotropic Fourier
// transform
double setFTQRIntJP(double r, double kp, double rmin1, double rmax1,int flag)
{
  static double k;
  static double rmin;
  static double rmax;
  double xi;
  double ans;

  if (flag == 1) {
    k = kp;
    rmin=rmin1;
    rmax=rmax1;
    return 0.0;
  }
  if(r>rmax) return 0.0;
  if(r<rmin) return 0.0;

  splint(rftJ,xiftJ,xiftspJ,nrJ,r,&xi);
  ans = xi/k;  //sin included in integration routine - r in spline
  //ans = xi/k*r;
  return ans;
}

// Dummy integrand for making power spectrum 
double FTQRIntJP(double r, void * params)
{
  return setFTQRIntJP(r,0.0,0.0,0.0,0);
}

///////////////////////////////////////////////////////////////////////////
// Calculate DD power spectraum
///////////////////////////////////////////////////////////////////////////
double Reionization::getCorrDD(double z, double rvec[], double xi[], int nr)
{
  int i;
  double r,corrDD;

  //calculate r values

  xidd(r,z,c,1);

  //calculate correlation function

  for(i=1;i<=nr;i++){
    r=rvec[i];
    corrDD=xidd(r,z,c,0);
    xi[i]=corrDD;
  }
}

double Reionization::getCorrXD(double z, double rvec[], double xi[], int nr)
{
  int i;
  double r,corrXD;
  double xion;

  xion=a->getXI(z);
  if((xion<XICUT)||(fabs(xion-1.0)<1.0e-3)){  // fully neutral or ionized
    for(i=1;i<=nr;i++) xi[i]=0.0;
  }else{
    //    xidd(r,z,c,1);
    for(i=1;i<=nr;i++){
      r=rvec[i];
      corrXD=correlateXD(r,z);
      xi[i]=corrXD;
      cout<<"XD \t"<<r<<"\t"<<corrXD<<endl;
    }
  }

  return 0.0;
}

double Reionization::getCorrXX(double z, double rvec[], double xi[], int nr)
{
  int i;
  double r,corrXX;
  double xion;

  //calculate correlation function
  xion=a->getXI(z);
  if((xion<XICUT)||(fabs(xion-1.0)<1.0e-3)){  // fully neutral or ionized
    for(i=1;i<=nr;i++) xi[i]=0.0;
  }else{
    //xidd(r,z,c,1);
    for(i=1;i<=nr;i++){
      r=rvec[i];
      corrXX=correlateXX(r,z);
      xi[i]=corrXX;
      cout<<"XX \t"<<r<<"\t"<<corrXX<<endl;
    }
  }
  
  return 0.0;
}

double Reionization::getCorrTB(double z, double rvec[], double xi[], int nr, double corrDD[], double corrXD[], double corrXX[], double xH)
{
  int i;
  double corr;


  if(xH<1.0e-3) return 0.0;

  for(i=1;i<=nr;i++){
    corr=corrXX[i]*(1.0+corrDD[i]);
    corr+= xH*xH*corrDD[i];
    corr+=corrXD[i]*(2.0*xH+corrXD[i]);
    xi[i]=corr/xH/xH;
  }

  return 0.0;
}


int Reionization::getRVec(double rvec[])
{
  int i,nr;
  int shift,shift2,shift3;
  double r,rmax,rmin,rstep;

  shift=20;
  shift2=40;
  shift3=20;

  nr=shift+shift2+shift3-1;

  //calculate r values

  rmin=1.0e-4;
  rmax=10.0;
  rstep=exp(log(rmax/rmin)/(double)(shift-1));
  r=rmin;
  for(i=1;i<=shift;i++){
    rvec[i]=r;
    r*=rstep;
  }
  
  rmin=r;
  rmax=300.0;
  rstep=exp(log(rmax/rmin)/(double)(shift2-1));
  for(i=shift+1;i<=shift+shift2;i++){
    rvec[i]=r;
    r*=rstep;
  }

  rmin=r;
  rmax=2200.0;
  rstep=exp(log(rmax/rmin)/(double)(shift3-1));
  for(i=shift+shift2+1;i<=nr;i++){
    rvec[i]=r;
    r*=rstep;
  }

  return nr;
}

///////////////////////////////////////////////////////////////////////
//  Mean clumping factor from Furlanetto and Oh (2005)
///////////////////////////////////////////////////////////////////////
//mean clumping factor including bubbles
//
double Reionization::clumpingBub(double z, double xi)
{
  double *ans;
  int goodSteps,badSteps;
  double step,oldans;
  double mMin,mMax;
  double clump;
  double lambda0,mLambda;

  setZeta(xi/c->fColl(z));

  ans = dvector(1,1);
  ans[1] = 0.0;
  setClumpingBubIntegrand(0.0,ans,ans,z,this,a,c,1);
  mMin=getZeta()*coolMass(c,z);  //minimum mass of bubble
  lambda0=6.0e6/(c->hubbleZ(z)*c->getH()*H0_CMSMPC)*(1.0+z);
  mLambda=4.0*PI/3.0*CRITDENMSOLMPC*c->getOm0hh()*pow(lambda0,3.0);
  cout<<mMin<<"\t"<<mLambda<<endl;
  //  mMin=fmax(mMin,mLambda);
  while (1) {
    step = mMin/10.0;
    mMax = mMin*10.0;
    oldans = ans[1];
    odeint(ans,1,mMin,mMax,1.0e-4,step,0.0,&goodSteps,&badSteps,
	   clumpingBubIntegrand,bsstep);
    if (fabs((ans[1]-oldans)/ans[1]) < 1.0e-4)
      break;
    mMin = mMax;
  }
  clump = ans[1]/(CRITDENMSOLMPC*c->getOm0hh());
  free_dvector(ans,1,1);

  return clump;
}

// Integrand for computing the bubble filling fraction.
//   m = bubble mass (Msun) [integration variable]
//   y,deriv = parameters used in the integration
//
void setClumpingBubIntegrand(double mb, double y[], double deriv[], 
			double z1, Reionization *rz1, Astrophysics *a1, Cosmology *c1, int flag)
{
  static double z;
  static Reionization *rz;
  static Astrophysics *a;
  static Cosmology *c;
  double weight;
  double deltai;
  double clumpBub,deltax;
  double *BB;
  double rb,lambda0,fv;
  double mfp_mean;

  if (flag == 1) {
    z = z1;
    rz = rz1;
    a=a1;
    c=c1;
    return;
  }
  BB=dvector(1,3);

  rb=pow(mb/(4.0*PI/3.0*CRITDENMSOLMPC*c->getOm0hh()),1.0/3.0);
  lambda0=6.0e6/(c->hubbleZ(z)*c->getH()*H0_CMSMPC)*(1.0+z);

  //if mfp at deltai=1.0 (mean density)
  //greater than bubble size take deltai=1.0
  fv=a->getVolumeMEHR(1.0,z);
  mfp_mean=lambda0/pow(1.0-fv,2.0/3.0);
  if(rb<mfp_mean){
    deltai=1.0;
  }else{
    fv=1.0-pow(lambda0/rb,1.5);
    deltai=a->deltaIMEHR(z,fv);
  }
  clumpBub=a->getClumpingMEHR(deltai,z);
 
  deltax=rz->barrierBub(mb,z,BB);
  //    deltax=0.0;

  weight=clumpBub*(1.0+deltax);
  deriv[1] = rz->nmBub(mb,z)*weight;
  //  cout<<mb<<"\t"<<deltai<<"\t"<<weight<<"\t"<<weight*rz->nmBub(mb,z)<<"\t"<<fv<<"\t"<<clumpBub<<endl;
  free_dvector(BB,1,3);
}

// Dummy function for integrating the bubble filling Q
void clumpingBubIntegrand(double m, double y[], double deriv[])
{
  setClumpingBubIntegrand(m,y,deriv,0.0,NULL,NULL,NULL,0);
}

// Find delta i using \lambda=R_b criterea
//
double Reionization::deltaIBub(double mb, double z)
{
  double deltai;
  double fv, lambda0;
  double rb;

  rb=pow(mb/(4.0*PI/3.0*CRITDENMSOLMPC*c->getOm0hh()),1.0/3.0)/(1.0+z); //physical bubble size
  lambda0=6.0e6/(c->hubbleZ(z)*c->getH()*H0_CMSMPC);
  
  fv=1.0-pow(lambda0/rb,1.5);
  //  cout<<fv<<"\t"<<rb<<"\t"<<lambda0<<endl;
  //be clever here to avoid writing more code
  //  if(fv<1.0e-4) return exp(-11.0);  //mfp>rb so can always grow
  deltai=a->deltaIMEHR(z,fv);
  return deltai;
}
