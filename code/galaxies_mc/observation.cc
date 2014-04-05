
#include <math.h>
#include <iostream>
#include <fstream>
#include "observation.h"
#include "dcosmology.h"
#include "dnumrecipes.h"
#include "astrophysics.h"
#include "twentyonecm.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//Constructor/Destructor
/////////////////////////////////////////////////////////////////////
Observation::Observation(Cosmology *c1)
{
  // cout <<"21cm constructor has been called." <<endl;
	c=c1;
	antennaeDistFlag=0;  //ugly hack to make sure reinitialise antennae dist
	foregroundFlag=1;
	optimiseFlag=0;
	fixed_area_flag=0;
  //cout << "21cm Constructor called successfully." <<endl; 
}

Observation::~Observation()
{
  // cout <<"Atomic destructor has been called." <<endl;
}

//////////////////////////////////////////////////////////////////////////
// Member Functions
/////////////////////////////////////////////////////////////////////////
void Observation::setCosmology(Cosmology *c1)
{
	c=c1;
}

double Observation::getDMin()
{
	return Dmin;
}

double Observation::getDMax()
{
	return Dmax;
}

double Observation::getAlpha()
{
	return alpha;
}

double Observation::getNAnt()
{
	return Nant;
}

double Observation::getNInnerRing()
{
	return nInnerRing;
}

double Observation::getRInnerRing()
{
	return rInnerRing;
}

double Observation::getFreqRes()
{
	return freqRes;
}

///////////////////////////////////////////////////////////////////////////
// Limits of Fisher integration
///////////////////////////////////////////////////////////////////////////
double Observation::getThetaMin(double z, double k)
{
	double y, cH;
	double zB, thetamin;

  y=surveyDepth(z);
  thetamin=acos(fmin(y*k/2.0/PI,1.0));

  //check for going into region where projected u_perp below minimum 
  //base lines have no sensitivity below this 
	double kmin ,x, u, r, thetab, lambda;
	x=surveyDistance(z);
	lambda=lambda21cm*(1.0+z);
	r=x*k*sin(thetamin)/2.0/PI*lambda;
	if(r<Dmin){
		thetab=asin(Dmin*2.0*PI/x/k/lambda);
		//cout<<thetamin<<"\t"<<thetab<<endl;
		thetamin=fmax(thetamin,thetab);
	}

	return thetamin;
}

double Observation::getThetaMax(double z, double k)
{
  double kstar,xres,x, cH;
  double thetamax;

  x=surveyDistance(z);
  xres=x*angularResolution(nu21cm/(1.0+z));
  kstar=2.0*PI/xres;  

  thetamax=asin(fmin(kstar/k,1.0));

	return thetamax;
}

double Observation::getSampleConstantD(double z)
{
	double D;
	double lambda,Aeff,cH,x,zB,y;

  lambda=lambda21cm*(1.0+z);
  Aeff=effectiveAntennaeArea(nu21cm/(1.0+z));
  x=surveyDistance(z);
  y=surveyDepth(z);

  //sample variance
  D=sqrt(4.0*PI*PI*Aeff/(x*x*y*lambda*lambda));


	return D/sqrt(Nfield);
}

double Observation::getThermalConstantE(double z)
{
	double E;
	double lambda,Aeff,cH,x,zB,y, Tsys;

  lambda=lambda21cm*(1.0+z);
  Tsys=tempSky(nu21cm/(1.0+z));
  Aeff=effectiveAntennaeArea(nu21cm/(1.0+z));
  x=surveyDistance(z);
  y=surveyDepth(z);

  //cosmic variance
  E=2.0*PI*sqrt(x*x*y)*lambda*lambda*lambda*Tsys*Tsys;
  E/=pow(Aeff,1.5)*B*Tint;

	return E/sqrt(Nfield);
}

//////////////////////////////////////////////////////////////////////////
// Initialise observations
/////////////////////////////////////////////////////////////////////////
void Observation::initObservation(double Nant1, double Atot1, double Dmin1, double Dmax1, double B1, double Tint1, double eta1,double alpha1, double Nfield1, double antennaeProp[])
{
	double fracInner, fracOuter1, fracOuter2;

	Nant=Nant1;
	Atot=Atot1;

	if(Atot<0.0){
		fixed_area_flag=1;
		Atot=fabs(Atot);
	}


	Dmin=Dmin1;
	Dmax=Dmax1;
	B=B1;
	Tint=Tint1;
	eta=eta1;
	alpha=alpha1;
	freqRes=8.0e3;
	Nfield=Nfield1;

	//set antennae properties
	fracInner=antennaeProp[1];
	fracOuter1=antennaeProp[2];
	fracOuter2=antennaeProp[3];

	nInnerRing=Nant*fracInner;
	nOuterRing1=Nant*fracOuter1;
	nOuterRing2=Nant*fracOuter2;

	rInnerRing=antennaeProp[4];
	rOuterRing1=antennaeProp[5];
	rOuterRing2=antennaeProp[6];

	Dmax=fmax(rInnerRing,fmax(rOuterRing1,rOuterRing2));
	Dmax*=2.0;

	ZOPTIMISE=8.0;   //typically optimse array for z=8.0;
                         //response drops at lower redshifts

}

///////////////////////////////////////////////////////////////////////////
// Init observation from file
///////////////////////////////////////////////////////////////////////////
void Observation::initObservationFromFile(string file)
{
	ifstream fin;
	string name;
	double Nant1, Atot1, Dmin1, Dmax1, B1, Tint1, eta1, alpha1, Nfield1;
	double *antennaeProp;

	antennaeProp=dvector(1,6);

	fin.open(file.c_str());
	fin>>Nant1;
	fin>>Atot1;
	fin>>Dmin1;
	fin>>Dmax1;
	fin>>B1;
	fin>>Tint1;
	fin>>eta1;
	fin>>alpha1;
	fin>>Nfield1;
	fin>>antennaeProp[1];
	fin>>antennaeProp[2];
	fin>>antennaeProp[3];
	fin>>antennaeProp[4];
	fin>>antennaeProp[5];
	fin>>antennaeProp[6];
	fin.close();

	Tint1*=60.0*60.0;  //convert from hours to seconds

	initObservation(Nant1,Atot1,Dmin1,Dmax1,B1,Tint1,eta1,alpha1,Nfield1,antennaeProp);

	free_dvector(antennaeProp,1,6);

}

void Observation::saveObservationToFile(string file)
{
	ofstream fout;

	fout.open(file.c_str());
	fout<<Nant;
	fout<<Atot;
	fout<<Dmin;
	fout<<Dmax;
	fout<<B;
	fout<<Tint/60.0/60.0;  //in hours
	fout<<eta;
	fout<<alpha;
	fout<<Nfield;
	fout<<nInnerRing/Nant;
	fout<<nOuterRing1/Nant;
	fout<<nOuterRing2/Nant;
	fout<<rInnerRing;
	fout<<rOuterRing1;
	fout<<rOuterRing2;
	fout.close();

}

///////////////////////////////////////////////////////////////////////////
// Observational performance
//////////////////////////////////////////////////////////////////////////

// Field of view
//
double Observation::fieldOfView(double nu)
{
	double lambda;
	lambda=SPEEDOFLIGHT_CGS/nu;
	return PI*pow(lambda/Dmin,2.0);
}

// Angular width of survey
//
double Observation::angularWidth(double nu)
{
	double lambda;
	lambda=SPEEDOFLIGHT_CGS/nu;
	return lambda/Dmin;
}

double Observation::surveyDepth(double z)
{
   double zP, zM, y, cH;
   cH=SPEEDOFLIGHT_CGS/c->getH()/UNH;
   zP=nu21cm/(nu21cm/(1.0+z)+B/2.0)-1.0;
   zM=nu21cm/(nu21cm/(1.0+z)-B/2.0)-1.0;
   y=fabs(c->confTime(zM)-c->confTime(zP))*cH/MPC;  //depth of volume probed

  return y;
}

double Observation::angularResolution(double nu)
{
	double lambda;
	lambda=SPEEDOFLIGHT_CGS/nu;
	return lambda/Dmax;
}

double Observation::surveyWidth(double z)
{
    double x,cH,theta;
    cH=SPEEDOFLIGHT_CGS/c->getH()/UNH;
    x=c->coordDistance(z)*cH/MPC;
    theta=sqrt(fieldOfView(nu21cm/(1.0+z))/PI);
    return x*theta; 
}

double Observation::surveyDistance(double z)
{
	double cH, x;
  cH=SPEEDOFLIGHT_CGS/c->getH()/UNH;
  x=c->coordDistance(z)*cH/MPC;

	return x;
}

double Observation::surveyVolume(double z)
{
	return pow(surveyWidth(z),2.0)*surveyDepth(z);
}

double Observation::getAntennaeArea()
{
	return Atot/Nant;
}

// Calculate the effective area of each antennae - this is frequency 
// dependent.
//   We assume that the physical collecting area is specified at z=8
// and scale from there to obtain the effective collecting area
//
double Observation::effectiveAntennaeArea(double nu)
{
  double Aeff,z;

  //fixed antennae size
  if(fixed_area_flag==1) return Atot/Nant;

  z=nu21cm/nu-1.0;
  Aeff=Atot/Nant*pow((1.0+z)/(1.0+ZOPTIMISE),2.0);  

  return Aeff;
}


// Foreground cutoff scale arising from finite frequency depth of the 
// observation.  Assumed that need frequency information to remove
// foregrounds.
//
// Units: kcut  :  Mpc^{-1}
//
double Observation::foregroundCutoff(double z)
{
	double cH, zB, y, kcut;

	y=surveyDepth(z);
        kcut=2.0*PI/y;

	return kcut;
}

//////////////////////////////////////////////////////////////////////////
// Antennae distribution
//////////////////////////////////////////////////////////////////////////

// Calculate antennae density
//
//  Units:  cm^{-2}
//
double Observation::antennaeDensity(double r)
{
	double rhomax;
	double rbreak;

	rhomax=1.0/Dmin/Dmin;

	if(antennaeDistFlag==0){
		//determine distribution of antennae
		

		cout<<"setting rbreak"<<endl;
		cout<<getAntennaeArea()<<"\t"<<Dmin*Dmin<<endl;
		rbreak=findDensityBreak();
		if(rbreak<0.0) rhomax=Nant/(PI*Dmax*Dmax/2.0/2.0);
		antennaeDistFlag=1;
		rbreakDist=rbreak;
	}

	rbreak=rbreakDist;
	if(r<rbreak){
		return rhomax;
	}else if(r>Dmax/2.0){
		return 0.0;
	}else{
		return rhomax*pow(rbreak/r,alpha);
	}

	cout<<"error in antennaeDensity"<<endl;
	return -1.0;
}

// Find break in antennae distribution to account for finite element size
//
double Observation::findDensityBreak()
{
	double tol(1.0e-4);
	double rbreak(Dmin);

	if(fabs(alpha)<1.0e-4) return -1.0;

	setDensityBreak(Dmin,this,1);
	rbreak=zriddrSimp(dummyDensityBreak,1.0e0,Dmax/2.0,tol);
	cout<<"rbreak="<<rbreak/100.0<<" m"<<endl;
	return rbreak;
}

double setDensityBreak(double rbreak, Observation *Obs1, int iflag)
{
	static Observation *Obs;
	double kernel;
	double rhomax,alpha;

	if(iflag==1){
		Obs=Obs1;
		return 0.0;
	}

	rhomax=1.0/Obs->getDMin()/Obs->getDMin();
	alpha=Obs->getAlpha();

	if(fabs(alpha-2.0)<1.0e-4){
		kernel=2.0*PI*rhomax*rbreak*rbreak*log(Obs->getDMax()/2.0/rbreak);	
	}else{
		kernel=2.0*PI*rhomax*pow(rbreak,alpha)/(2.0-alpha);
		kernel*=pow(Obs->getDMax()/2.0,2.0-alpha)-pow(rbreak,2.0-alpha);
	}
	kernel+=PI*rhomax*rbreak*rbreak;

	return Obs->getNAnt()-kernel;
}

double dummyDensityBreak(double r)
{
	return setDensityBreak(r,NULL,0);	
}

//code to calculate antennae distribution with possibility of 
// two outer rings with constant antennae densities

double Observation::antennaeDensityFull(double r)
{
	double rho, rhomax;
	double rbreak;

	rhomax=1.0/Dmin/Dmin;

	if(antennaeDistFlag==0){
		cout<<"setting rbreak"<<endl;
		rbreak=findDensityBreakFull();
		cout<<"rbreak="<<rbreak/100.0<<" m"<<endl;
		antennaeDistFlag=1;
		rbreakDist=rbreak;
	}

	rbreak=rbreakDist;
	if(rbreak<0.0){
		 rhomax=nInnerRing/(PI*rInnerRing*rInnerRing);
		 rbreak=rInnerRing;
	}


	if(r<rbreak){
		return rhomax;
	}else if(r<rInnerRing){
		return rhomax*pow(rbreak/r,alpha);	
	}else if(r<rOuterRing1){
		rho=nOuterRing1/PI/(rOuterRing1*rOuterRing1-rInnerRing*rInnerRing);
		return rho;
	}else if(r<rOuterRing2){
		rho=nOuterRing2/PI/(rOuterRing2*rOuterRing2-rOuterRing1*rOuterRing1);
		return rho;
	}else{
		return 0.0;
	}

	cout<<"error in antennaeDensity"<<endl;
	return -1.0;

}


// Find break in antennae distribution to account for finite element size
//
double Observation::findDensityBreakFull()
{
	double tol(1.0e-4);
	double rbreak(Dmin);

	if(fabs(alpha)<1.0e-4) return -1.0;  //flat distribution

	setDensityBreakFull(Dmin,this,1);
	rbreak=zriddrSimp(dummyDensityBreakFull,1.0e0,rInnerRing,tol);
	return rbreak;
}

double setDensityBreakFull(double rbreak, Observation *Obs1, int iflag)
{
	static Observation *Obs;
	double kernel;
	double rhomax,alpha;

	if(iflag==1){
		Obs=Obs1;
		return 0.0;
	}

	rhomax=1.0/Obs->getDMin()/Obs->getDMin();
	alpha=Obs->getAlpha();

	if(fabs(alpha-2.0)<1.0e-4){
		kernel=2.0*PI*rhomax*rbreak*rbreak*log(Obs->getRInnerRing()/rbreak);	
	}else{
		kernel=2.0*PI*rhomax*pow(rbreak,alpha)/(2.0-alpha);
		kernel*=pow(Obs->getRInnerRing(),2.0-alpha)-pow(rbreak,2.0-alpha);
	}
	kernel+=PI*rhomax*rbreak*rbreak;

	return Obs->getNInnerRing()-kernel;
}

double dummyDensityBreakFull(double r)
{
	return setDensityBreakFull(r,NULL,0);	
}

//Calculate baseline number density in (u,v,eta) space
//
double Observation::baselineDensityUV(double u, double z)
{
	double r,nuv;
	double lambda;

	lambda=lambda21cm*(1.0+z);
	r=u*lambda;
	nuv=getXCol(r);         // units cm^{-2}
	nuv*=lambda*lambda;     //dimensionless number density

	return nuv;
}

// Optimise antennae distribution for maximum sensitivity at redshift z
// Essentially need to set minimum baseline to lambda/2
//
void Observation::optimiseAntennae(double z)
{
	double lambda, Aeff;
	double ndipole;
	double ncore;
	int concentrate_flag(0);

	if(fabs(z-ZOPTIMISE)<1.0e-4){
		cout<<"antennae already optimised"<<z<<"\t"<<ZOPTIMISE<<endl;
		return;
	}
	
	cout<<"optimising antennae"<<endl;
	cout<<"old Dmin="<<Dmin<<endl;
	cout<<"old Aeff="<<effectiveAntennaeArea(nu21cm/(1.0+z))<<endl;

	//calculate number of dipoles per tile - basic things that fixes experiment
	ndipole=Atot/Nant;
	ndipole/=pow(lambda21cm*(1.0+ZOPTIMISE)/2.0,2.0);
	cout<<"n_dipole="<<ndipole<<endl;

	lambda=lambda21cm*(1.0+z);
	Dmin=sqrt(ndipole)*lambda/2.0;                 //set A_tile=A_eff
	Aeff=Dmin*Dmin;
	Atot=Nant*Aeff;                  //make sure total area matches
	ZOPTIMISE=z;                     //remember optimised for this z
	antennaeDistFlag=0;              //ensure recalculate antennae dist

	ncore=PI*pow(rInnerRing/Dmin,2.0); //max number of close packed antennae
	if(ncore<nInnerRing){
		cout<<"enlarging array to fit in all antennae"<<endl;
		cout<<"old Dmax="<<Dmax<<endl;
		//enforce max baseline twice that needed to fit in antennae
		rInnerRing=Dmin*sqrt(nInnerRing/PI)*2.0;
		Dmax=2.0*rInnerRing;
		cout<<"new Dmax="<<Dmax<<endl;
	}


	if(concentrate_flag==1 && ZOPTIMISE>3.0){
	//concentrate antennae so that filling fraction=1 throughout
	cout<<"old Dmax="<<Dmax<<endl;
	double trash;
	alpha=0.0;
	Dmax=2.0*sqrt(1.1*Atot/PI);
	rInnerRing=Dmax/2.0;
	cout<<"new Dmax="<<Dmax<<endl;
	}

	cout<<"new Dmin="<<Dmin<<endl;
	cout<<"new Aeff="<<Aeff<<endl;
}


//////////////////////////////////////////////////////////////////////////
//  Visibility distribution
//////////////////////////////////////////////////////////////////////////

//Calculate cross-correllation between antennae elements
//
//  Units:  n  :  cm^{-2}     (number density of visibilities)
//          b  :  cm          (baseline)
//
double Observation::xcolAntennae(double b)
{
	double n;
	double tol(1.0e-5);
	double rmin,rmax;

	if(b<Dmin) return 0.0;
	if(b>Dmax) return 0.0;  // no correlation beyond 2.0* (Dmax/2.0)

	rmin=Dmin;       //not sure on this lower limit
	rmax=Dmax/2.0;

	if(b>Dmax/2.0) rmin=b-Dmax/2.0;  //need b-r<Dmax/2 to get signal

	setXColAntennae(b,0.0,this,1);
	setRXCol(b,this,1);

	if(Dmax>2.0*rInnerRing*1.01){
		n=qromb(dummyRXCol,Dmin,rInnerRing,tol);	
		n+=qromb(dummyRXCol,rInnerRing,rOuterRing1,tol);
		n+=qromb(dummyRXCol,rOuterRing1,rmax,tol);
	}else{
		n=qromb(dummyRXCol,rmin,rmax,tol);
	}

	return n;
}

double setXColAntennae(double rp, double phi, Observation *Obs1, int iflag)
{
	static Observation *Obs;
	static double b;
	double xcol,d;	

	if(iflag==1){
		Obs=Obs1;
		b=rp;
		return 0.0;
	}

	double phicut;
	phicut=((pow(Obs->getDMax()/2.0,2.0)-b*b-rp*rp)/2.0/b/rp);
	d=sqrt(b*b+rp*rp+2.0*b*rp*cos(phi));
	xcol=Obs->antennaeDensityFull(rp)*Obs->antennaeDensityFull(d)*rp;
//	if(xcol==0) cout<<b<<"\t"<<rp<<"\t"<<Obs->getDMax()/2.0<<"\t"<<phicut<<"\t"<<acos(phicut)<<"\t"<<phi<<"\t"<<d<<"\t"<<xcol<<endl;

	return xcol;
}

double rpsaveXC;

double dummyPhiXCol(double phi)
{
	return setXColAntennae(rpsaveXC,phi,NULL,0);
}

double dummyRXCol(double r)
{
//	double tol(1.0e-5);
//	rpsaveXC=r;
//	return qromb(dummyPhiXCol,-PI,PI,tol);
	return setRXCol(r,NULL,0);
}

double setRXCol(double r, Observation *Obs1, int iflag)
{
	static double b;
	static Observation *Obs;
	double phimin, phimax;
	double tol(1.0e-5);
	double phicut;

	if(iflag==1){
		b=r;
		Obs=Obs1;
		return 0.0;
	}

	phicut=((pow(Obs->getDMax()/2.0,2.0)-b*b-r*r)/2.0/b/r);

	rpsaveXC=r;
	phimin=0.0;
	phimax=2.0*PI;
	if(fabs(phicut)<1.0){
		 phicut=1.001*acos(phicut);
		 phimin=phicut;
		 phimax=2.0*PI-phicut;
	}

//	cout<<PI<<"\t"<<phicut<<"\t"<<phimax<<endl;
	return qromb(dummyPhiXCol,phimin,phimax,tol);
}

//Calculation of xcol is slow, so do it once and then use a spline
void Observation::initXColSpline()
{
	double r,rstep;
	int i,nstep(50);
	nstep=75;
	int nlow,nhigh;

	double *rvec, *xcvec;

	rvec=dvector(1,nstep);
	xcvec=dvector(1,nstep);

//	if(Dmax>2.0*rInnerRing){
	if(Dmax>1.0e30){
		nlow=25;
		nhigh=nstep-nlow;

		rstep=exp(log(rInnerRing/Dmin)/(double)(nlow-1));
		r=Dmin;
		for(i=1;i<=nlow;i++){
			rvec[i]=log(r);
			xcvec[i]=xcolAntennae(r);
			cout<<Dmax<<"\t"<<exp(rvec[i])<<"\t"<<xcvec[i]<<endl;
			r*=rstep;
		}

		rstep=exp(log(Dmax/r)/(double)(nhigh-1));
		for(i=nlow+1;i<=nstep;i++){
			rvec[i]=log(r);
			xcvec[i]=xcolAntennae(r);
			//cout<<Dmax<<"\t"<<exp(rvec[i])<<"\t"<<xcvec[i]<<endl;
			r*=rstep;
		}

	}else{
		rstep=exp(log(Dmax/Dmin)/(double)(nstep-1));
		r=Dmin;
		for(i=1;i<=nstep;i++){
			rvec[i]=log(r);
			xcvec[i]=xcolAntennae(r);
			//cout<<Dmax<<"\t"<<exp(rvec[i])<<"\t"<<xcvec[i]<<endl;
			r*=rstep;
		}
	}

	visDens.setSplineSP(nstep,rvec,xcvec);

	free_dvector(rvec,1,nstep);
	free_dvector(xcvec,1,nstep);
}

double Observation::getXCol(double r)
{
	if(!visDens.checkInit() || antennaeDistFlag==0){
		cout<<"intialising xcol spline"<<endl;
		initXColSpline();
//		forceNormVisibilities();
		cout<<"initialisation complete"<<endl;
	}

//	if(r<exp(visDens.getXMin())) return 0.0;
//	if(r>exp(visDens.getXMax())) return 0.0;

	if(r<exp(visDens.getXMin())) return visDens.returnValue(visDens.getXMin());
	if(r>exp(visDens.getXMax())) return visDens.returnValue(visDens.getXMax());

	return visDens.returnValue(log(r));
}

//check normalization of visibilities
//integral should give N_vis=Nant*(Nant-1)
double Observation::numberVisibilities()
{
	double nvisInt;
	double tol(1.0e-4);

	setNumberVisibilitiesInt(0.0,this,1);
	nvisInt=PI*qromb(dummyNumberVisibilitiesInt,Dmin,Dmax,tol);

	return nvisInt;
}

double setNumberVisibilitiesInt(double r, Observation *Obs1, int iflag)
{
	static Observation *Obs;
	double visdens;

	if(iflag==1){
		Obs=Obs1;
		return 0.0;
	}

	visdens=Obs->getXCol(r)*r;

	return visdens;
}

double dummyNumberVisibilitiesInt(double r)
{
	return setNumberVisibilitiesInt(r,NULL,0);
}

void Observation::forceNormVisibilities()
{
	double ntrue, ncalc;
	xColRenorm=1.0;

	ncalc=numberVisibilities();
	ntrue=Nant*(Nant-1.0)/2.0;

	xColRenorm=ntrue/ncalc;
	cout<<"renorm: "<<xColRenorm<<endl;

	if(fabs(xColRenorm-1.0)>1.0) cout<<"Warning: xColRenorm large!"<<endl;
}

/////////////////////////////////////////////////////////////////////////
// Antennae response function
/////////////////////////////////////////////////////////////////////////

// Individual dipoles have a response function in real space 
// that looks like a dipole (unsurprisingly)
//
//Units:   theta  : radians
//          nu	  : Hz
double Observation::antennaAngularResponse(double theta, double nu)
{
	double w;
	double fwhm;

	fwhm=angularWidth(nu);
	cout<<fwhm*180/PI<<endl;
	if(theta<fwhm){
		w=pow(cos(PI*theta/2.0/fwhm),2.0);
	}else{
		w=0.0;
	}
	return w;
}



///////////////////////////////////////////////////////////////////////////
// 21 cm observations
//////////////////////////////////////////////////////////////////////////

//sky temperature
// Units: Tsky Kelvin
double Observation::tempSky(double nu)
{
  double nu0(1.8e8);
  double T0(180.0);
  double Tsky;

  Tsky=T0*pow(nu/nu0,-2.6);

// for testing against mcquinn calculation
  double z;
  z=nu21cm/nu-1.0;

  if(MCQUINN==1){
  if(fabs(z-8.0)<1.0e-3){
	return 440.0;
  }else if(fabs(z-6.0)<1.0e-3){
	return 250.0;
  }else if(fabs(z-12.0)<1.0e-3){
	return 1000.0;
  }
	}

  
  return Tsky;
}

//Calculate the error on a measurement of the power spectrum
// 
// Assume log steps \Delta k = eta * k
//
double Observation::powerSpectrumError(double z, double k, double ps)
{
  double Nbase;
  double D, E, x, y, lambda, Tsys, Aeff;
  double kstar;
  double cH;
  double errP;
  double zB;
 double xstar,yfreq;
  double thetamin,thetamax,tol(1.0e-4);
  double errInt;
  double xres;

  if(optimiseFlag==1) optimiseAntennae(z);

  if(k<foregroundCutoff(z) && foregroundFlag==1) {
	cout<<"erased by foreground"<<endl;
	return 1.0e15;
  }

  Nbase=Nant*(Nant-1.0)/2.0;
  lambda=lambda21cm*(1.0+z);
  Tsys=tempSky(nu21cm/(1.0+z));

  Aeff=effectiveAntennaeArea(nu21cm/(1.0+z));

  cH=SPEEDOFLIGHT_CGS/c->getH()/UNH;
  zB=nu21cm/(nu21cm/(1.0+z)+B)-1.0;

  x=surveyDistance(z);
  y=surveyDepth(z);

  E=2.0*PI*sqrt(x*x*y)*lambda*lambda*lambda*Tsys*Tsys;
  E/=sqrt(eta)*pow(Aeff,1.5)*B*Tint;
//  E*=2.0;  //different convention between Bowman and McQuinn

  D=sqrt(4.0*PI*PI*Aeff/(x*x*y*lambda*lambda*eta));

  xres=x*angularResolution(nu21cm/(1.0+z));

  kstar=2.0*PI/xres;  

  setPowerSpectrumErrorKernel(0.0,k,z,x,D,E,ps,this,1);
  thetamin=acos(fmin(y*k/2.0/PI,1.0));
  thetamax=asin(fmin(kstar/k,1.0));
//  cout<<k<<"\t"<<thetamin<<"\t"<<thetamax<<"\t"<<kstar<<"\t"<<2.0*PI/y<<"\t"<<kstar/k<<"\t"<<y*k/2.0/PI<<endl;
  errInt=qromb(dummyPowerSpectrumErrorKernel,thetamin,thetamax,tol);
  errP=1.0/sqrt(k*k*k*fabs(errInt));
  errP/=sqrt(Nfield);

  if(errP>1.0e30) return 1.0e30;  //handle infiinites
  return errP;
}


// 
double setPowerSpectrumErrorKernel(double theta, double k1, double z1, double x1, double D1, double E1, double ps1, Observation *Obs1, int iflag)
{
	static double x,k, D, E, ps, z;
	static Observation *Obs;
	double kernel;
	double n;

	if(iflag==1){
		k=k1;
		z=z1;
		x=x1;
		D=D1;
		E=E1;
		ps=ps1;
		Obs=Obs1;
		return 0.0;
	}

	n=Obs->baselineDensityUV(x*k*sin(theta)/2.0/PI,z);
	
//	cout<<x<<"\t"<<k<<"\t"<<n<<"\t"<<z<<"\t"<<theta<<"\t"<<x*k*sin(theta)/2.0/PI<<endl;
	kernel=sin(theta);
	kernel*=pow(n/(D*ps*n+E),2.0);
	return kernel;
}

double dummyPowerSpectrumErrorKernel(double theta)
{
	return setPowerSpectrumErrorKernel(theta,0.0,0.0,0.0,0.0,0.0,0.0,NULL,0);
}

//Calculate the error on a measurement of the angular components of the
// power spectrum P_{mu^0}, P_{mu^2}, and P_{mu^3}
// 
// Assume log steps \Delta k = eta * k
//
double Observation::powerSpectrumErrorMu(double z, double k, double ps0, double ps2, double ps4, int order)
{
	double **fisher, **ifisher;
	double **identity;
	double errP;
	int index;
	int i,j;

	index=1+order/2;

	if(optimiseFlag==1) optimiseAntennae(z);

	fisher=dmatrix(1,3,1,3);
	ifisher=dmatrix(1,3,1,3);
	identity=dmatrix(1,3,1,3);

	if(k<foregroundCutoff(z)&& foregroundFlag==1) {
		cout<<"erased by foreground"<<endl;
		return 1.0e15;
  	}

	//initialise fisher matrix
	fisher[1][1]=fisherErrorMu(z,k,ps0,ps2,ps4,0);
	fisher[1][2]=fisherErrorMu(z,k,ps0,ps2,ps4,2);
	fisher[1][3]=fisherErrorMu(z,k,ps0,ps2,ps4,4);
	fisher[2][1]=fisher[1][2];
	fisher[2][2]=fisher[1][3];
	fisher[2][3]=fisherErrorMu(z,k,ps0,ps2,ps4,6);
	fisher[3][1]=fisher[1][3];
	fisher[3][2]=fisher[2][3];
	fisher[3][3]=fisherErrorMu(z,k,ps0,ps2,ps4,8);

	//check to see if matrix is singular
	if(fabs(detMatrix(fisher,3))<1.0e-40){
	cout<<"singular fisher matrix"<<endl;
	return 1.0e30;	
	}

	//calculate inverse fisher matrix
	invertMatrix(fisher,3,ifisher);

	multMatrix(fisher,ifisher,identity,3);

	double idmax;
	idmax=0.0;
//	cout<<"identity"<<endl;
	for(i=1;i<=3;i++){
		for(j=1;j<=3;j++){
			if(i==j) identity[i][j]-=1.0;
//			cout<<identity[i][j]<<"\t";
			if(fabs(identity[i][j])>idmax) idmax=fabs(identity[i][j]);
		}
//		cout<<endl;
	}
	cout<<endl;
	if(idmax>1.0e-6) cout<<"Fisher problems idmax="<<idmax<<endl;	
	

	if(ifisher[index][index]<0.0){
		 cout<<"diagonal element of ifisher negative"<<endl;
	}

	errP=sqrt(fabs(ifisher[index][index]));
 	errP/=sqrt(Nfield);

	free_dmatrix(fisher,1,3,1,3);
	free_dmatrix(ifisher,1,3,1,3);
	free_dmatrix(identity,1,3,1,3);

	return errP;
}

// Calculate elemenets of angular dependent Fisher matrix for calculating
// the the errors on angular dependent components of power spectrum
//
double Observation::fisherErrorMu(double z, double k, double ps0, double ps2, double ps4, int order)
{
  double Nbase;
  double D, E, x, y, lambda, Tsys, Aeff;
  double kstar;
  double cH;
  double errP;
  double zB;
 double xstar,yfreq;
  double thetamin,thetamax,tol(1.0e-4);
  double errInt;
  double xres;

  Nbase=Nant*(Nant-1.0)/2.0;
  lambda=lambda21cm*(1.0+z);
  Tsys=tempSky(nu21cm/(1.0+z));
  Aeff=effectiveAntennaeArea(nu21cm/(1.0+z));

  cH=SPEEDOFLIGHT_CGS/c->getH()/UNH;

  x=surveyDistance(z);
  y=surveyDepth(z);

  //cosmic variance
  E=2.0*PI*sqrt(x*x*y)*lambda*lambda*lambda*Tsys*Tsys;
  E/=sqrt(eta)*pow(Aeff,1.5)*B*Tint;

  //thermal noise
  D=sqrt(4.0*PI*PI*Aeff/(x*x*y*lambda*lambda*eta));

  xres=x*angularResolution(nu21cm/(1.0+z));
  zB=nu21cm/(nu21cm/(1.0+z)+freqRes)-1.0;
  yfreq=fabs(c->confTime(z)-c->confTime(zB))*cH/MPC;

  kstar=2.0*PI/xres;  

  setFisherErrorMuKernel(0.0,k,z,x,D,E,ps0,ps2,ps4,order,this,1);
  thetamin=acos(fmin(y*k/2.0/PI,1.0));
  thetamax=asin(fmin(kstar/k,1.0));

  errInt=qromb(dummyFisherErrorMuKernel,thetamin,thetamax,tol);
  errP=k*k*k*errInt;

  return errP;
}


// 
double setFisherErrorMuKernel(double theta, double k1, double z1, double x1, double D1, double E1, double ps01, double ps21, double ps41, int order1, Observation *Obs1, int iflag)
{
	static double x,k, D, E, ps0, ps2, ps4, z;
	static int order;
	static Observation *Obs;
	double kernel;
	double ps, mu;
	double n;

	if(iflag==1){
		k=k1;
		z=z1;
		x=x1;
		D=D1;
		E=E1;
		ps0=ps01;
		ps2=ps21;
		ps4=ps41;
		order=order1;
		Obs=Obs1;
		return 0.0;
	}

	n=Obs->baselineDensityUV(x*k*sin(theta)/2.0/PI,z);
	mu=cos(theta);	

	ps=ps0+ps2*pow(mu,2.0)+ps4*pow(mu,4.0);

//	cout<<x<<"\t"<<k<<"\t"<<n<<"\t"<<z<<"\t"<<theta<<"\t"<<x*k*sin(theta)/2.0/PI<<endl;
	kernel=sin(theta);
	kernel*=pow(mu,(double)(order));
	kernel*=pow(n/(D*ps*n+E),2.0);
	return kernel;
}

double dummyFisherErrorMuKernel(double theta)
{
	return setFisherErrorMuKernel(theta,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,NULL,0);
}

///////////////////////////////////////////////////////////////////
//return inverse of matrix leaving original untouched
//
// Uses gaussj - somewhat unreliable
//
void invertMatrix2(double **fisher,int nparam, double **ifisher)
{
  int i,j,m(1);
  double **b;
  
  b=dmatrix(1,nparam,1,m);

  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      ifisher[i][j]=fisher[i][j];
    }
  }

  for(i=1;i<=nparam;i++){
    b[i][1]=1.0;
  }
  //need this routine from NR
  gaussj(ifisher,nparam,b,m);

  free_dmatrix(b,1,nparam,1,m);
}

//save a matrix to file
void writeMatrix(double **M,int nparam, char *file)
{
  ofstream fout;
  int i,j;

  fout.open(file);
  fout.precision(NPRECISION);
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      if(j==nparam){
	fout<<M[i][j]<<endl;	
      }else{
	fout<<M[i][j]<<"\t";
      }
    }
  }

  fout.close();
} 

void writeMatrixS(double **M,int nparam, string file)
{
  ofstream fout;
  int i,j;

  fout.open(file.c_str());
  fout.precision(NPRECISION);
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      if(j==nparam){
	fout<<M[i][j]<<endl;	
      }else{
	fout<<M[i][j]<<"\t";
      }
    }
  }

  fout.close();
} 

double detMatrix(double **M, int nparam)
{
	int i,j, *indx;
	double d;

	if(nparam>3){
		cout<<"fine if non-singular"<<endl;
		ludcmp(M,nparam,indx,&d);	
		for(j=1;j<=nparam;j++) d*=M[j][j];
	}else if(nparam==2){
		//2 by 2 matrix
		d=M[1][1]*M[2][2]-M[1][2]*M[2][1];
	}else{
		//3 by 3 matrix
		d=M[1][1]*(M[2][2]*M[3][3]-M[3][2]*M[2][3]);
		d+= -M[1][2]*(M[2][1]*M[3][3]-M[2][3]*M[3][1]);
		d+=M[1][3]*(M[2][1]*M[3][2]-M[2][2]*M[3][1]);
	}

	return d;
}

// Multiply the square N*N matrices A and B to obtain C
//
void multMatrix(double **A, double **B, double **C, int N)
{
	int i,j,k;

	for(i=1;i<=N;i++){
		for(j=1;j<=N;j++){
			C[i][j]=0.0;
			for(k=1;k<=N;k++){
				C[i][j]+=A[i][k]*B[k][j];
			}
		}
	}

}

//return inverse of matrix leaving original untouched
//
// Uses LU decomposition - fairly reliable even for nearly singular matrices
//
void invertMatrix(double **fisher,int nparam, double **ifisher)
{
  int i,j;
  double **a, d, *col;
  int *indx;

  a=dmatrix(1,nparam,1,nparam);
  col=dvector(1,nparam);
  indx=ivector(1,nparam);

  //write fisher into a since inversion destroys original matrix
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      a[i][j]=fisher[i][j];
    }
  }

  //perform inversion

  ludcmp(a,nparam,indx,&d);

  for(j=1;j<=nparam;j++){
	for(i=1;i<=nparam;i++) col[i]=0.0;
	col[j]=1.0;
	lubksb(a,nparam,indx,col);
	for(i=1;i<=nparam;i++) ifisher[i][j]=col[i];
  }

  free_dmatrix(a,1,nparam,1,nparam);
  free_dvector(col,1,nparam);
  free_ivector(indx,1,nparam);
}

//return inverse of matrix leaving original untouched
//
// Uses SVD decomposition - note clear that its an improvement on LU decomp
//
void invertMatrixSVD(double **fisher,int nparam, double **ifisher)
{
  int i,j;
  double wmax, wmin, **a, **u, *w, **v, *x, *col;

  a=dmatrix(1,nparam,1,nparam);
  u=dmatrix(1,nparam,1,nparam);
  v=dmatrix(1,nparam,1,nparam);
  x=dvector(1,nparam);
  w=dvector(1,nparam);
  col=dvector(1,nparam);

  //write fisher into a since inversion destroys original matrix
  for(i=1;i<=nparam;i++){
    for(j=1;j<=nparam;j++){
      a[i][j]=fisher[i][j];
    }
  }

  //perform inversion

  svdcmp(a,nparam,nparam,w,v);
  wmax=0.0;
  for(j=1;j<=nparam;j++) if(w[j]>wmax) wmax=w[j];
  //set threshold
  wmin=1.0e-12*wmax;
  for(j=1;j<=nparam;j++){
	if(w[j]<wmin){
		cout<<"zeroing\t"<<j<<"\t"<<w[j]<<endl;
		w[j]=0.0;
	}
  }

  // calculate inverse matrix
  for(j=1;j<=nparam;j++){
	for(i=1;i<=nparam;i++) col[i]=0.0;
	col[j]=1.0;
	svbksb(a,w,v,nparam,nparam,col,x);
	for(i=1;i<=nparam;i++) ifisher[i][j]=x[i];
  }

  free_dmatrix(a,1,nparam,1,nparam);
  free_dmatrix(u,1,nparam,1,nparam);
  free_dmatrix(v,1,nparam,1,nparam);
  free_dvector(x,1,nparam);
  free_dvector(col,1,nparam);
  free_dvector(w,1,nparam);
}


//return inverse of matrix leaving original untouched
// USES ZERO OFFSET
// Uses LU decomposition - fairly reliable even for nearly singular matrices
//
void invertMatrixZO(double **fisher,int nparam, double **ifisher)
{
  int i,j;
  double **a, d, *col;
  int *indx;

  a=dmatrix(1,nparam,1,nparam);
  col=dvector(1,nparam);
  indx=ivector(1,nparam);

  //write fisher into a since inversion destroys original matrix
  for(i=0;i<=nparam-1;i++){
    for(j=0;j<=nparam-1;j++){
      a[i+1][j+1]=fisher[i][j];
    }
  }

  //perform inversion

  ludcmp(a,nparam,indx,&d);

  for(j=1;j<=nparam;j++){
	for(i=1;i<=nparam;i++) col[i]=0.0;
	col[j]=1.0;
	lubksb(a,nparam,indx,col);
	for(i=1;i<=nparam;i++) ifisher[i-1][j-1]=col[i];
  }

  free_dmatrix(a,1,nparam,1,nparam);
  free_dvector(col,1,nparam);
  free_ivector(indx,1,nparam);
}

////////////////////////////////////////////////////////////////////////
// Numerical recipes gaussj routine ///////////////////////////////////
///////////////////////////////////////////////////////////////////////

#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI

/////////////////////////////////////////////////////////////////
