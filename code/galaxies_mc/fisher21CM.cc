/* fisherCMB.cc
 * Contains functions for calculating fisher matrix for a Galaxy experiment
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include "fisher21CM.h"
#include "fisher.h"
#include "dcosmology.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include "dnumrecipes.h"
#include "spline.h"

using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 Fisher21CM::Fisher21CM()
{
	//cout<<"Fisher21CM constructor called"<<endl;
	fisher_flag=0;
       
}

 Fisher21CM::~Fisher21CM()
{
	//cout<<"Fisher21CM destructor called"<<endl;

}

//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////

//initialise Galaxy survey parameters
void Fisher21CM::initObservation(string name1, double z1, double kmin1, double kmax1, Observation *Obs1, int regime_flag)
{
	exp_name=name1;
	exp_z=z1;

	kmin=kmin1;
	kmax=kmax1;

//	dirbase="./data/cosmology";
	dirbase=dirbase+"_z"+numToString(exp_z);

	SurveyObs=Obs1;

	post_reion_flag=regime_flag;   //=1 for post reionization regime

	setSurveyConstants();

	cout<<"survey: "<<exp_name<<" initialised"<<endl;
}

//////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////

string Fisher21CM::getSurveyName()
{
	return exp_name;
}

double Fisher21CM::getSurveyZ()
{
	return survey_z;
}

double Fisher21CM::getKMin()
{
	return fmax(SurveyObs->foregroundCutoff(exp_z)+0.001,kmin);
}

double Fisher21CM::getKMax()
{
	return kmax;
}

void Fisher21CM::setKMax(double kmax1)
{
	kmax=kmax1;
}

////////////////////////////////////////////////////////////////////////////
// Bias parameters
/////////////////////////////////////////////////////////////////////////////

// Calculate amplitude scaling for 21 cm power spectrum -- crudely this
// is a scale independent bias after reionization and T_b before
//
double Fisher21CM::getBias(CosmParam *c)
{
	double tb(0.015);

	if(post_reion_flag) return 0.03;
	
	tb=26.0*(c->ombhh/0.022);   //in mK
	tb*=sqrt(0.15/c->ommhh*(1.0+exp_z)/10.0);
	return tb;
}

////////////////////////////////////////////////////////////////////////////
// Set useful survey constants to avoid frequent recalculation
/////////////////////////////////////////////////////////////////////////////
void Fisher21CM::setSurveyConstants()
{
	xSurveyVolume=SurveyObs->surveyDistance(exp_z);
	sampleConstantD=SurveyObs->getSampleConstantD(exp_z);
	thermalConstantE=SurveyObs->getThermalConstantE(exp_z);
}

void Fisher21CM::setNoThermalNoise()
{
	thermalConstantE=0.0;
}


////////////////////////////////////////////////////////////////////////////
//  Set galaxy survey parameter files
////////////////////////////////////////////////////////////////////////////
double Fisher21CM::setGalaxyFiles()
{
  int i;
  //parameters for calculating derivatives
  for(i=0;i<=nparamFC;i++) setGalaxyPair(i);
  return 1.0;
}


double Fisher21CM::setGalaxyPair(int iflag)
{
  double p_step;
  CosmParam c_p;
  CosmParam c_m;
  string base,name,file,tagp,tagm;
  string tag;
  double vary(0.05);
  double bias,ascal2;
  int i;
  int noncosm_flag(0);  //indicates noncosmological parameter

  base=dirbase;
  tagp="_param_p.dat";
  tagm="_param_m.dat";

  	tag=param_tags[iflag];
	name=tag;

	//for non-cosmological parameters use fiducial
	if(name=="_bias") noncosm_flag=1;
	if(name=="_lda") noncosm_flag=1;
	if(name=="_lh") noncosm_flag=1;
	if(name=="_lg") noncosm_flag=1;
	if(name=="_ascal2") noncosm_flag=1;  //handles better this way

	if(noncosm_flag==1){
		name="_fiducial";
	}

	file=base+name+tagp;
	c_p=loadCosmParam(file);
	bias=fiducialFC.bias;
	c_p.bias=bias;

	file=base+name+tagm;
	c_m=loadCosmParam(file);
	c_m.bias=bias;

	//handle distance parameters
	if(distance_flag==0){
		c_m.lda=log(angDiamDistance(exp_z,&c_m));
		c_p.lda=log(angDiamDistance(exp_z,&c_p));
		c_p.lh=log(hubbleZ(exp_z,&c_p));
		c_m.lh=log(hubbleZ(exp_z,&c_m));
		c_m.lg=log(growthZ(exp_z,&c_m));
		c_p.lg=log(growthZ(exp_z,&c_p));
        }else if(distance_flag==1 && tag=="_fiducial"){
		c_m.lda=log(angDiamDistance(exp_z,&c_m));
		c_p.lda=log(angDiamDistance(exp_z,&c_p));
		c_p.lh=log(hubbleZ(exp_z,&c_p));
		c_m.lh=log(hubbleZ(exp_z,&c_m));
		c_m.lg=log(growthZ(exp_z,&c_m));
		c_p.lg=log(growthZ(exp_z,&c_p));

	}else{
		c_m.lda=fiducialFC.lda;
		c_p.lda=fiducialFC.lda;	
		c_m.lh=fiducialFC.lh;
		c_p.lh=fiducialFC.lh;
		c_m.lg=fiducialFC.lg;
		c_p.lg=fiducialFC.lg;
	}
		

	//handle derivatives with respect to non-cosmological parameters
	if(tag=="_fiducial"){
		bias=getBias(&c_p);
		c_p.bias=bias;
		c_m.bias=bias;
		fiducialFC=c_p;
	}else if(tag=="_bias"){
		param_values[iflag]=bias;  //insert fiducial value of bias
		p_step=vary*bias;
		c_p.bias+=p_step;
		c_m.bias-=p_step;
	}else if(tag=="_lda"){
		param_values[iflag]=c_p.lda;  //insert fiducial value of lda
		p_step=vary*c_p.lda;
		c_p.lda+=p_step;
		c_m.lda-=p_step;
	}else if(tag=="_lg"){
		param_values[iflag]=c_p.lg;  //insert fiducial value of lg
		p_step=vary*c_p.lg;
		c_p.lg+=p_step;
		c_m.lg-=p_step;
	}else if(tag=="_lh"){
		param_values[iflag]=c_p.lh;  //insert fiducial value of lh
		p_step=vary*c_p.lh;
		c_p.lh+=p_step;
		c_m.lh-=p_step;
	}else if(tag=="_ascal2"){
		param_values[iflag]=c_p.ascal2;  //insert fiducial value of lh
		p_step=vary*c_p.ascal2;
		c_p.ascal2+=p_step;
		c_m.ascal2-=p_step;
	}


	//handle file related details for non-cosmological parameters
	if(noncosm_flag==1){
		c_p.file=base+tag+"_p.dat";
		c_m.file=base+tag+"_m.dat";

		c_p.name=base+tag+"_param_p.dat";
		c_m.name=base+tag+"_param_m.dat";

		//replicate fiducial transfer and matter spectrum - wasteful
		name="cp "+dirbase+"_fiducial_p.dat"+" "+dirbase+tag+"_p.dat";
		system(name.c_str()); 
		name="cp "+dirbase+"_fiducial_m.dat"+" "+dirbase+tag+"_m.dat";
		system(name.c_str());

		base=base+"_transfer";
		name="cp "+base+"_fiducial_p.dat"+" "+base+tag+"_p.dat";
		system(name.c_str()); 
		name="cp "+base+"_fiducial_m.dat"+" "+base+tag+"_m.dat";
		system(name.c_str()); 
	}


  saveCosmParam(&c_m,c_m.name);
  saveCosmParam(&c_p,c_p.name);

  return 1.0;
}

/////////////////////////////////////////////////////////////////////////
// Calculate Galaxy Fisher Matrix
/////////////////////////////////////////////////////////////////////////

//global variables for fisher matrix integration
double lkint21CM;
//double lkminint,lkmaxint;

//Evaluate the Fisher Matrix via the Tegmark approximation
// Generates the fisher matrix in fisher[nmat][nmat] given
// the survey volume (vsurvey), the galaxy number density (density)
// and a fiducial cosmology (c)
void Fisher21CM::fisherMatrix21CM()
{
  double tol(1.0e-5);
  double fM;
  double lkminint,lkmaxint;
  string file_dlPi;
  string file_dlPj;

  int i,j;
  ofstream fout;
  ofstream fsave;
	string file_save;
  string file;
  string file_fid;
  string tag;
  int nmat;
  double **fisher, **ifisher;
  double *ans;
  int goodSteps,badSteps;
  double step;
	double lkmin, lkmax, lkstep;
	int b;

	ans=dvector(1,1);

  nparamFC=param_tags.size()-1;
  nmat=nparamFC;
//  for(i=0;i<=nparamFC;i++) cout<<param_tags[i]<<endl;
  string file_deriv[nmat+1];

  fisher=dmatrix(1,nmat,1,nmat);
  ifisher=dmatrix(1,nmat,1,nmat);

  cout<<"setting galaxy files"<<endl;
  setGalaxyFiles();

  file_fid=param_tags[0];

  for(i=0;i<=nparamFC-1;i++) file_deriv[i]=param_tags[i+1];

  cout<<"calculating fisher matrix"<<endl;
  //everything is done using files to create the fisher matrix
  //I'll exploit this to automate the proceedure

    lkminint=log(getKMin());
    lkmaxint=log(getKMax());

  for(i=0;i<nmat;i++){
    for(j=0;j<=i;j++){
      file_dlPi=file_deriv[i];
      file_dlPj=file_deriv[j];

//	file_save="./fisherint"+file_dlPi+file_dlPj+".dat";
//	fsave.open(file_save.c_str());	

       setFisherKernelK21CM(0.0,NULL,NULL,exp_z,SurveyObs,1);
       setFisherKernel21CM(0.0,NULL,NULL,file_fid,this,file_dlPi,file_dlPj,0);
	lkstep=0.1;
	ans[1]=0.0;
	lkmin=lkminint;
	lkmax=lkmin+lkstep;
	b=0;
	while(1){
	 	if(lkmax>lkmaxint){
			b=1;
			lkmax=lkmaxint;
		}
	step=0.01*(lkmax-lkmin);
  odeint(ans,1,lkmin,lkmax,tol,step,0.0,&goodSteps,&badSteps,getFisherKernelK21CM,bsstep);
	lkmin=lkmax;
	lkmax+=lkstep;
//	fsave<<exp(lkmin)<<"\t"<<exp(lkmax)<<"\t"<<ans[1]<<endl;
	if(b==1) break;
	}
	fM=ans[1];
      setFisherKernel21CM(0.0,NULL,NULL,file_fid,this,file_dlPi,file_dlPj,2);
      fisher[i+1][j+1]=fM;
      if(i!=j) fisher[j+1][i+1]=fM;
   cout <<file_deriv[i]<<"\t"<<file_deriv[j]<<"\t"<<fM<<endl;
//	fsave.close();
    }
  }

  //save information about fisher matrix and inverse into class
  setFisherMatrix(fisher,nmat,&fiducialFC);

  //write result to file
  tag="_z"+numToString(exp_z)+".dat";
  saveFisher(tag);

  //output error information
  cout <<"Galaxy survey errors: "<<exp_name<<endl;

  for(i=1;i<=nparamFC;i++){
	cout<<"delta"<<param_tags[i]<<"\t"<<param_values[i]<<"\t"<<sqrt(ifisherFC[i][i])<<"\t"<<sqrt(ifisherFC[i][i])/param_values[i]<<endl;
  }

  free_dmatrix(fisher,1,nmat,1,nmat);
  free_dmatrix(ifisher,1,nmat,1,nmat);
  free_dvector(ans,1,1);

}

void getFisherKernelK21CM(double lk, double y[], double deriv[])
{
  setFisherKernelK21CM(lk,y,deriv,0.0,NULL,0);
}

void setFisherKernelK21CM(double lk, double y[], double deriv[], double z1, Observation *Obs1, int iflag)
{
  double gFKM;
  double tol(1.0e-5);
  static Observation *Obs;
  static double z;
  double thetamin, thetamax;
  double thetaminint, thetamaxint;
  double thetastep;
  double *ans;
  int goodSteps,badSteps;
  double step;

  int b;

  if(iflag==1){
	z=z1;
	Obs=Obs1;
	return;
  }
	ans=dvector(1,1);

  lkint21CM=lk;

  thetaminint=Obs->getThetaMin(z,exp(lk));
  thetamaxint=Obs->getThetaMax(z,exp(lk));

//  gFKM=qromb(getFisherKernel21CM,thetamin,thetamax,tol);

	thetamin=thetaminint;
	thetastep=(thetamaxint-thetaminint)/4.9; //avoid tmax to tmax integral
	thetamax=thetamin+thetastep;
	b=0;
	ans[1]=0.0;
  while(1){
	if(thetamax>thetamaxint){
		thetamax=thetamaxint;
		b=1;
	}

	  step=(thetamax-thetamin)/40.0;
  odeint(ans,1,thetamin,thetamax,tol,step,0.0,&goodSteps,&badSteps,getFisherKernel21CM,bsstep);
//	cout<<exp(lk)<<"\t"<<thetamin<<"\t"<<thetamax<<"\t"<<thetaminint<<"\t"<<thetamaxint<<"\t"<<ans[1]<<"\t"<<lk<<endl;
	thetamin=thetamax;
	thetamax+=thetastep;
	if(b==1) break;
  }
	gFKM=ans[1];
	deriv[1]=gFKM;
  //	return ans[1];
  //return gFKM;	
}


void getFisherKernel21CM(double theta, double y[], double deriv[])
{
  string blank;
  setFisherKernel21CM(theta,y,deriv,blank,NULL,blank,blank,1);
}

void setFisherKernel21CM(double theta, double y[], double deriv[], string file_fid, Fisher21CM *survey1,string file_dlPi1in, string file_dlPj1in,int iflag)
{
  static CosmParam c;
  static Fisher21CM *survey;
  double sFK;
  double k,lk;
  static CosmParam ci_m, ci_p;
  static CosmParam cj_m, cj_p;
  double dlPi,dlPj;
  string file,base,tagm,tagp;
  static string file_dlPi1, file_dlPj1;

  tagm="_param_m.dat";
  tagp="_param_p.dat";


  lk=lkint21CM;
  k=exp(lk);

  if(iflag==0){  //initiate
    survey=survey1;
    base=survey->returnDirbase();

    file_dlPi1=file_dlPi1in;
    file_dlPj1=file_dlPj1in;


    file=base+file_fid+tagm;
    c=survey->loadCosmParam(file);
    file=base+file_dlPi1+tagm;
    ci_m=survey->loadCosmParam(file);
    file=base+file_dlPi1+tagp;
    ci_p=survey->loadCosmParam(file);
    file=base+file_dlPj1+tagm;
    cj_m=survey->loadCosmParam(file);
    file=base+file_dlPj1+tagp;
    cj_p=survey->loadCosmParam(file);

	ci_m.daratio=survey1->daratio(&ci_m);
	ci_p.daratio=survey1->daratio(&ci_p);
	ci_m.hratio=survey1->hratio(&ci_m);
	ci_p.hratio=survey1->hratio(&ci_p);
	cj_m.daratio=survey1->daratio(&cj_m);
	cj_p.daratio=survey1->daratio(&cj_p);
	cj_m.hratio=survey1->hratio(&cj_m);
	cj_p.hratio=survey1->hratio(&cj_p);

    survey->setDerivSplines(file_dlPi1,file_dlPj1,1);
    return;
  }

  if(iflag==2){  //clean up spline routines
	    survey->setDerivSplines(file_dlPi1,file_dlPj1,2);
    return;
  }

  //calculate power spectrum derivatives
  dlPi=survey->derivePower21CM(k,theta,ci_p,ci_m,file_dlPi1,1); 
  dlPj=survey->derivePower21CM(k,theta,cj_p,cj_m,file_dlPj1,2);

  //if difference is very small the above method becomes numerically unstable
  if((dlPi!=dlPi)||(dlPj!=dlPj)) cout<<k<<"\t"<<"Error in setFisherKernel"<<endl;

  sFK=dlPi*dlPj*survey->variance21CM(k,theta);
  sFK*=k*k*k;
  sFK*=sin(theta);  //angular weighting
  deriv[1]=sFK;
}

//Calculate variance on measurement for fiducial cosmology
double Fisher21CM::variance21CM(double k, double theta)
{
  double E,D,n,ps,var,z,x;
  
  z=exp_z;
  x=xSurveyVolume;

  n=SurveyObs->baselineDensityUV(x*k*sin(theta)/2.0/PI,z);
  ps=getPowerSpectrumObs(k,theta,&fiducialFC,&splineFid);
  E=thermalConstantE;
  D=sampleConstantD;

 if(E<1.0e-30){
	var=pow(1.0/D/ps,2.0);
 }else{
	 var=pow(n/(D*ps*n+E),2.0);
 }

//  cout<<k<<"\t"<<theta<<"\t"<<z<<"\t"<<E<<"\t"<<D<<"\t"<<n<<"\t"<<ps<<"\t"<<var<<endl;

  return var;
}

//Calculate effective volume of survey for fiducial cosmology
double Fisher21CM::effectiveVolume(double k, double theta)
{
  double ps,var;
  
  ps=getPowerSpectrumObs(k,theta,&fiducialFC,&splineFid);

 var=variance21CM(k,theta);
 var*=2.0*pow(2.0*PI,2.0);
 var*=ps*ps;

//  cout<<k<<"\t"<<theta<<"\t"<<z<<"\t"<<E<<"\t"<<D<<"\t"<<n<<"\t"<<ps<<"\t"<<var<<endl;

  return var;
}


///////////////////////////////////////////////////////////////////////
//  Code to calculate the galaxy power spectrum
///////////////////////////////////////////////////////////////////////

//calculate the derivative of power sepctrum
double Fisher21CM::derivePower21CM(double k, double theta, CosmParam c_p, CosmParam c_m, string tag, int iflag)
{
	double dP, p_step;
	Spline *splineP, *splineM;

	if(iflag==1){
		splineP=&splinePow1p;
		splineM=&splinePow1m;
	}else if(iflag==2){
		splineP=&splinePow2p;
		splineM=&splinePow2m;
	}	

    p_step=getParamFromTag(tag,&c_p)-getParamFromTag(tag,&c_m);

  dP=((getPowerSpectrumObs(k,theta,&c_p,splineP))-(getPowerSpectrumObs(k,theta,&c_m,splineM)))/(p_step);

	return dP;
} 

//set splines for power spectrum derivatives
void Fisher21CM::setDerivSplines(string tag1, string tag2, int iflag)
{
 	string base,tagm,tagp;
	string transfertag;

	tagm="_m.dat";
  	tagp="_p.dat";
	transfertag="_transfer";

	//set splines
	if(iflag==1){
	    base=dirbase;

	assignSpline(&splineFid,base,"_fiducial"+tagp);
	assignSpline(&splinePow1p,base,tag1+tagp);
	assignSpline(&splinePow1m,base,tag1+tagm);
	assignSpline(&splinePow2p,base,tag2+tagp);
	assignSpline(&splinePow2m,base,tag2+tagm);

	return;
	}

	//clean splines
	if(iflag==2){
	splineFid.cleanSpline();	
	splinePow1p.cleanSpline();
	splinePow1m.cleanSpline();
	splinePow2p.cleanSpline();
	splinePow2m.cleanSpline();

	return;
	}

	cout<<"Must enter 1 or 2 for iflag"<<endl;
	return;

}


//Calculate the galaxy power spectrum for the parameters specified in c
//at the observed values of (kref,muref) using a spline transfer function 
double Fisher21CM::getPowerSpectrumObs(double kref, double theta,CosmParam *c, Spline *spline)
{
  double ps;
  double pshot;
  double mu,muref;

  double alcockp;
  double kpara2_ref,kperp2_ref;
  double kpara2,kperp2;
  double k,kref2;
  double hratio, daratio;

  pshot=0.0;
  daratio=exp(c->lda)/exp(fiducialFC.lda);
  hratio=exp(c->lh)/exp(fiducialFC.lh);

  muref=cos(theta);

  //relate reference k,mu to the true k,mu
  kref2=kref*kref;
  kpara2_ref=kref2*muref*muref;
  kperp2_ref=kref2*(1.0-muref*muref);
  kperp2=kperp2_ref/pow(daratio,2.0);
  kpara2=kpara2_ref*pow(hratio,2.0);
  k=sqrt(kperp2+kpara2);
  mu=sqrt(kpara2/kpara2+kperp2);

  alcockp=hratio/pow(daratio,2.0);  //Alcock-Paczynski factor

  ps=powerSpectrum(k,c,spline);
  ps*=pow(c->bias/1.0e3,2.0); 		      //bias
  ps*=1.0+2.0*pow(mu,2.0)+pow(mu,4.0);        //peculiar velocities
  ps*=alcockp;

  return ps+pshot;
}


double Fisher21CM::powerSpectrumError(double k)
{
	double tol(1.0e-5);
	double thetamin, thetamax;
	double errInt;
	double z;
	string base,file;

	z=exp_z;

	if(k<SurveyObs->foregroundCutoff(z)) return 1.0e30;

  thetamin=SurveyObs->getThetaMin(z,k);
  thetamax=SurveyObs->getThetaMax(z,k);

    file=dirbase+"_fiducial_param_p.dat";;
    fiducialFC=loadCosmParam(file);

	fiducialFC.lda=log(angDiamDistance(exp_z,&fiducialFC));
	fiducialFC.lh=log(hubbleZ(exp_z,&fiducialFC));

	    base=dirbase; 
	assignSpline(&splineFid,base,"_fiducial_p.dat");

	setKErrorKernel21CM(k,0.0,this,1);
	errInt=qromb(getKErrorKernel21CM,thetamin,thetamax,tol);
	return 1.0/sqrt(fabs(errInt));
}

double setKErrorKernel21CM(double k1, double theta, Fisher21CM *survey1, int iflag)
{
	static double k;
	static Fisher21CM *survey;
	double sFK;
	double eta(0.5); //k step size dk=eta*k

	if(iflag==1){
		k=k1;
		survey=survey1;
		return 0.0;
	}

  sFK=survey->variance21CM(k,theta);
  sFK*=k*k*k*eta;
  sFK*=sin(theta);  //angular weighting

	return sFK;
}

double getKErrorKernel21CM(double theta)
{
	return setKErrorKernel21CM(0.0,theta,NULL,0);
}


void Fisher21CM::outPowerErrors()
{
	double bias;
	double k, ps, err, veff;
	string file;
	ofstream fout;

	file=exp_name+"_k_errors.dat";

	bias=fiducialFC.bias/1.0e3;

	fout.open(file.c_str());
	k=kmin;
	while(k<kmax){
		ps=powerSpectrumFid(k,exp_z)*bias*bias;
		err=powerSpectrumError(k);
		veff=averagedEffectiveVolume(k);
		fout<<k<<"\t"<<ps<<"\t"<<err<<"\t"<<bias<<"\t"<<veff<<endl;
		k*=1.2;
	}
	fout.close();
	
}


double Fisher21CM::averagedEffectiveVolume(double k)
{
	double tol(1.0e-5);
	double thetamin, thetamax;
	double errInt;
	double z;
	string base,file;

	z=exp_z;

	if(k<SurveyObs->foregroundCutoff(z)) return 1.0e-30;

  thetamin=SurveyObs->getThetaMin(z,k);
  thetamax=SurveyObs->getThetaMax(z,k);

    file=dirbase+"_fiducial_param_p.dat";;
    fiducialFC=loadCosmParam(file);

	fiducialFC.lda=log(angDiamDistance(exp_z,&fiducialFC));
	fiducialFC.lh=log(hubbleZ(exp_z,&fiducialFC));

	    base=dirbase;
	assignSpline(&splineFid,base,"_fiducial_p.dat");

	setKVolumeKernel21CM(k,0.0,this,1);
	errInt=qromb(getKVolumeKernel21CM,thetamin,thetamax,tol);
	return errInt;
}

double setKVolumeKernel21CM(double k1, double theta, Fisher21CM *survey1, int iflag)
{
	static double k;
	static Fisher21CM *survey;
	double sFK;
	double eta(0.5); //k step size dk=eta*k

	if(iflag==1){
		k=k1;
		survey=survey1;
		return 0.0;
	}

  sFK=survey->effectiveVolume(k,theta);
  sFK*=sin(theta);  //angular weighting
  sFK/=2.0;  //averaging over mu

	return sFK;
}

double getKVolumeKernel21CM(double theta)
{
	return setKVolumeKernel21CM(0.0,theta,NULL,0);
}