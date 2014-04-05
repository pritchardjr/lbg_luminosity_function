/* fisherLF.cc
 * Contains functions for calculating fisher matrix for a Galaxy experiment
  */

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "fisherLF.h"
#include "galaxies.h"
#include "haloDensity.h"


using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

 FisherLF::FisherLF()
{
	//cout<<"FisherLF constructor called"<<endl;
  	fisher_flag=0;     

        min_mag= -22.0;  //minimum absolute magnitude reached v.bright
}

 FisherLF::~FisherLF()
{
	//cout<<"FisherLF destructor called"<<endl;

}

//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////

//Units
//volume Mpc^{-3}
//maglim - absolute magnitude limit
void FisherLF::initExperiment(string name, double zmin, double zmax, double anglex1, double angley1, double maglim1, int nfield1, int nbin1)
{
   double area;
	exp_name=name;
        exp_z=(zmin+zmax)/2.0;
        min_z=zmin;
        max_z=zmax;
        anglex=anglex1;
        angley=angley1;
        max_mag=maglim1;
        nfield=nfield1;
        nbin=nbin1;

        area=anglex*angley/60.0/60.0;
        volume=volumeFromArea(zmin,zmax,area)*nfield;

        cout<<"area="<<area*60.0*60.0<<"\t"<<"volume="<<volume<<endl;

        magbin=0.25;
        nbin=ceil((max_mag-min_mag)/magbin);
        min_mag=max_mag-nbin*magbin;

	cout<<"survey: "<<exp_name<<" initialised"<<endl;
}

//////////////////////////////////////////////////////////////////////
// Member functions
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Routines for accessing non-cosmology parameters in nice way
//////////////////////////////////////////////////////////////////////
// LF has many non-cosmological parameters and want to access
// them in a natural way while still using the overall CosmParam structure

void FisherLF::setLFParam(double phi0, double M0, double alpha0)
{
  vector<double> reion_param;
  vector<string> reion_tags;
  string name;

  reion_param.clear();
  reion_param.push_back(phi0);
  reion_param.push_back(M0);
  reion_param.push_back(alpha0);

  reion_tags.clear();
  reion_tags.push_back("_lphi0");
  reion_tags.push_back("_M0");
  reion_tags.push_back("_alpha0");

  setSupParams(reion_param);
  setSupParamsTags(reion_tags);

  addParam("_lphi0");
  addParam("_M0");
  addParam("_alpha0");

  //simplistic way of getting the fiducial file set up
  name=dirbase+"_fiducial"+"_param_p.dat";
  saveCosmParam(&fiducialFC,name);
  name=dirbase+"_fiducial"+"_param_m.dat";
  saveCosmParam(&fiducialFC,name);

}
////////////////////////////////////////////////////////////////////////////
//  Set galaxy survey parameter files
////////////////////////////////////////////////////////////////////////////
double FisherLF::setGlobalFiles()
{
  int i;
  //parameters for calculating derivatives
  for(i=0;i<=nparamFC;i++) setGlobalPair(i);
  return 1.0;
}


double FisherLF::setGlobalPair(int iflag)
{
  double p_step(0.0);
  CosmParam c_p;
  CosmParam c_m;
  string base,name,file,tagp,tagm,tag;
  double vary(0.05);
  //vary=0.01;
  int noncosm_flag(0);  //indicates noncosmological parameter

  base=dirbase;
  tagp="_param_p.dat";
  tagm="_param_m.dat";

  	tag=param_tags[iflag];
	name=tag;
        //cout<<"Calculating history for: "<<tag<<endl;

	//for non-cosmological parameters use fiducial
        //if(tag.find("_sparam",0)!=string::npos) noncosm_flag=1;

        if(tag.find("_phi0",0)!=string::npos) noncosm_flag=1;
        if(tag.find("_M0",0)!=string::npos) noncosm_flag=1;
        if(tag.find("_alpha0",0)!=string::npos) noncosm_flag=1;        
        
	if(noncosm_flag==1){
		name="_fiducial";
	}

//hack
        name="_fiducial";

	file=base+name+tagp;
	//cout<<file<<endl;
	c_p=loadCosmParam(file);
        if(fiducialFC.sparam.size()>0){
           c_p=incSupParams(fiducialFC.sparam,c_p);
           c_p=incSupParamsTags(fiducialFC.sparam_tags,c_p);
        }

	file=base+name+tagm;
	//cout<<file<<endl;
	c_m=loadCosmParam(file);
        if(fiducialFC.sparam.size()>0){
           c_m=incSupParams(fiducialFC.sparam,c_m);
           c_m=incSupParamsTags(fiducialFC.sparam_tags,c_m);
        }

	//handle derivatives with respect to non-cosmological parameters
	if(tag=="_fiducial"){
		fiducialFC=c_p;
	}else{
           string tags;
           tags=tag;
           tags.erase(0,7);  
           param_values[iflag]=getParamFromTag(tag,&c_p);
           p_step=vary*getParamFromTag(tag,&c_p)/2.0;
	   if(fabs(p_step)<1.0e-6) p_step=1.0e-4;  //avoid zero
           //c_p.sparam[atoi(tags.c_str())]+=p_step;
           //c_m.sparam[atoi(tags.c_str())]-=p_step;

           //cout<<"B: "<<tag<<"\t"<<getParamFromTag(tag,&c_p)<<"\t"<<p_step<<endl;
           setParamFromTag(getParamFromTag(tag,&c_p)+p_step,tag,&c_p);
           //cout<<"A: "<<tag<<"\t"<<getParamFromTag(tag,&c_p)<<endl;

           //cout<<"B: "<<tag<<"\t"<<getParamFromTag(tag,&c_m)<<"\t"<<p_step<<endl;
           setParamFromTag(getParamFromTag(tag,&c_m)-p_step,tag,&c_m);
           //cout<<"A: "<<tag<<"\t"<<getParamFromTag(tag,&c_m)<<endl;
        }

	//handle file related details for non-cosmological parameters
        //if(reset_flag==1){
		c_p.file=base+tag+"_p.dat";
		c_m.file=base+tag+"_m.dat";

		c_p.name=base+tag+"_param_p.dat";
		c_m.name=base+tag+"_param_m.dat";

		c_p.cl_file=base+tag+"_cl_p.dat";
		c_m.cl_file=base+tag+"_cl_m.dat";

		c_p.filebase=tag;
		c_m.filebase=tag;
                c_p.p_step=p_step;
                c_m.p_step=p_step;
                //}
        
        if(noncosm_flag==1 && model_flag!=0){
		//replicate fiducial transfer and matter spectrum - wasteful
		name="cp "+dirbase+"_fiducial_p.dat"+" "+dirbase+tag+"_p.dat";
		system(name.c_str()); 
		name="cp "+dirbase+"_fiducial_m.dat"+" "+dirbase+tag+"_m.dat";
		system(name.c_str()); 
	}


        //cout<<"E: "<<tag<<"\t"<<getParamFromTag(tag,&c_m)<<endl;
        //  cout<<"E: "<<tag<<"\t"<<getParamFromTag(tag,&c_p)<<endl;
  saveCosmParam(&c_m,c_m.name);
  saveCosmParam(&c_p,c_p.name);

  return 1.0;
}


//////////////////////////////////////////////////////////////////////
// Initialisation routine
//////////////////////////////////////////////////////////////////////
void FisherLF::fisherMatrix()
{
  double fM, fMij;
  string file_dlPi;
  string file_dlPj;

  long i, j,l, m;
  ofstream fout;
  ofstream fsave;
  string file_save;
  string file_fid;
  int nmat;
  double **fisher, **ifisher;

  double magmin, magmax;
  double *maggrid;
  double magstep;
  int npoints(nbin);  //number grid points in each dimension
  //convergent to 1e-5 level for npoints=1000
  double **variancegrid;
  double **derivativegrid;
  CosmParam c_p, c_m, c;
  string file, base, tagm, tagp, tag;

  tagm="_m.dat";
  tagp="_p.dat";
  base=returnDirbase();

  cout<<nbin<<endl;

  setGlobalFiles();

  //assign memory for arrays
  nparamFC=param_tags.size()-1;
  nmat=nparamFC;
  string file_deriv[nmat+1];
  fisher=dmatrix(1,nmat,1,nmat);
  ifisher=dmatrix(1,nmat,1,nmat);

  maggrid=dvector(0,npoints-1);
  variancegrid=dmatrix(0,npoints-1,0,npoints-1);
  derivativegrid=dmatrix(0,npoints-1,1,nparamFC);
  //derivSgrid=d3tensor(0,npoints-1,0,npoints-1,1,nparamFC);

  //set up of files needed for fisher matrix calculation

  file_fid=param_tags[0];

  for(i=0;i<=nparamFC-1;i++) file_deriv[i]=param_tags[i+1];

  cout<<"calculating fisher matrix"<<endl;

  //begin calculation of fisher matrix

  //determine grid in uperp-kpara space for summation
  // grid has log spacing but put in actual value to array
  magmin=min_mag;
  magmax=max_mag;

  cout<<"Limits: "<<magmin<<"\t"<<magmax<<endl;

  magstep=(magmax-magmin)/(npoints-1);

  for(i=0;i<npoints;i++)  maggrid[i]=magmin+magstep*i;


  //set up fiducial spline
  file=base+"_fiducial"+"_param"+tagm;
  c=loadCosmParam(file);
  
  //if(model_flag==1) assignSpline(&splineFid,base,"_fiducial"+tagm);
  if(model_flag==1) splineFid.loadFileSpline(base+"_fiducial"+tagm,6,1,2); //astro
  if(model_flag==2) splineFid.loadFileSpline(base+"_fiducial"+tagm,3,1,3); //turn
  cout<<base<<endl;


  //calculate variance grid in uv space
  cout<<"Calculating variance grid"<<endl;
  for(i=0;i<npoints;i++){
     for(j=0;j<npoints;j++){
        cout<<i<<"\t"<<j<<endl;
        variancegrid[i][j]=getVariance(maggrid[i],maggrid[j]);
     }
  }

  //now begin calculation of derivative grids
  cout<<"calculating derivative grids"<<endl;
  for(l=1;l<=nparamFC;l++){
    //set up parameters and splines for derivative
    tag=file_deriv[l-1];
    cout<<tag<<endl;
    file=base+tag+"_param"+tagm;
    c_m=loadCosmParam(file);
    file=base+tag+"_param"+tagp;
    c_p=loadCosmParam(file);
    
    if(model_flag==1) splinePow1p.loadFileSpline(base+tag+tagp,6,1,2); //astro
    if(model_flag==1) splinePow1m.loadFileSpline(base+tag+tagm,6,1,2); //astro
    if(model_flag==2) splinePow1p.loadFileSpline(base+tag+tagp,3,1,3); //turn
    if(model_flag==2) splinePow1m.loadFileSpline(base+tag+tagm,3,1,3); //turn

    //loop over uv space
    for(i=0;i<npoints;i++){
       derivativegrid[i][l]=deriveNbin(maggrid[i],c_p,c_m,tag,1);
    }
      if(model_flag==1) splinePow1p.cleanSpline();
      if(model_flag==1) splinePow1m.cleanSpline();
  }

  //want inverse of Dij for fisher matrix calculation
  invertMatrixZO(variancegrid,npoints,variancegrid);

  //now actually calculate the fisher matrix by direct summation
  for(l=1;l<=nparamFC;l++){
    for(m=1;m<=l;m++){
      fM=0.0;
      //summation loops
      for(i=0;i<npoints;i++){
         for(j=0;j<npoints;j++){
            //first term from LF derivatives
            fMij=derivativegrid[i][l]*derivativegrid[j][m]*variancegrid[i][j];
            //cout<<l<<"\t"<<m<<"\t"<<i<<"\t"<<j<<"\t"<<derivativegrid[i][l]<<"\t"<<derivativegrid[i][m]<<"\t"<<variancegrid[i][j]<<endl;;
            //cout<<l<<"\t"<<m<<"\t"<<variancegrid[i]<<endl;;
            //cout<<l<<"\t"<<derivativegrid[i][m]<<endl;
            fM+=fMij;

            //second term from sample variance derivatives
            fMij=0.0;

            fM+=fMij;
         }      
      } 

      fisher[l][m]=fM;
      if(l!=m) fisher[m][l]=fM;
      //cout<<file_deriv[l-1]<<"\t"<<file_deriv[m-1]<<"\t"<<fM<<endl;
    }
  }

  //done with calculation now handle data in clever ways
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
  free_dvector(maggrid,0,npoints-1);
  free_dmatrix(variancegrid,0,npoints-1,0,npoints-1);
  free_dmatrix(derivativegrid,0,npoints-1,1,nparamFC);
}

double FisherLF::getVariance(double mi, double mj)
{
  double var;
  double Sij, Pij;

  //Sample variance contribution to covariance
  Sij=getSampleVariance(mi,mj,&fiducialFC);

  //Poisson contribution to the covariance
  if(fabs(mi-mj)<magbin/2.0){
     Pij=Nbin(mi,&fiducialFC)/volume;
  }else{
     Pij=0.0;
  }


  var=Sij+Pij;

  //cout<<mi<<"\t"<<mj<<"\t"<<var<<endl;
  return var;
}

double FisherLF::getSampleVariance(double mi, double mj, CosmParam *c)
{
  double Sij;
  double ni, nj, bi, bj;
  double window;

  //return 0.0;
  //although this runs the bias calculation is horribly slow.
  //since for fixed cosmology we're just looking up the value of the bias
  //and number of halos as a function of mass many times should spline these
  //two quantities. Loses generality of being able to take derivatives wrt 
  // cosmology though - or at least only lose speed advantage if done right

  ni=Nbin(mi,&fiducialFC);
  nj=Nbin(mj,&fiducialFC);
  cout<<"N: "<<ni<<"\t"<<nj<<endl;
  //cout<<mi<<"\t"<<magbin<<"\t"<<exp_z<<endl;
  bi=biasM(mi,magbin,exp_z,&fiducialFC);
  bj=biasM(mj,magbin,exp_z,&fiducialFC);
  cout<<"B: "<<bi<<"\t"<<bj<<endl;
  window=getWindow(c);
  cout<<"W="<<window<<endl;

  Sij=ni*nj*bi*bj/nfield*window;

  return Sij;
}

double FisherLF::deriveS(double mi, double mj,CosmParam c_p,CosmParam c_m,string tag,int iflag)
{
  double deriv, p_step;
  p_step=getParamFromTag(tag,&c_p)-getParamFromTag(tag,&c_m);

  deriv=(getSampleVariance(mi,mj,&c_p)-getSampleVariance(mi,mj,&c_m))/p_step;

  return deriv;
}

double FisherLF::deriveNbin(double m,CosmParam c_p,CosmParam c_m,string tag,int iflag)
{
  double deriv, p_step;
  p_step=getParamFromTag(tag,&c_p)-getParamFromTag(tag,&c_m);

  deriv=(Nbin(m,&c_p)-Nbin(m,&c_m))/p_step;

  return deriv;
}

//number of galaxies in a given magnitude bin
double FisherLF::Nbin(double M, CosmParam *c, Spline *spline)
{
   double lphi0,phi0, M0, alpha0;
   double nmean;
   Galaxies Gal;
   int i,npoints(10);
   double Mmin, Mmax, Mstep;

  //get parameters from tags specified when reioniztion model is set
  lphi0=getParamFromTag("_lphi0",c);
  M0=getParamFromTag("_M0",c);
  alpha0=getParamFromTag("_alpha0",c);
  phi0=exp(lphi0*log(10.0));

  //integrate LF over magnitude bin to get number
  //NEED TO LOOK AT THIS INTEGRATION IN MORE DETAIL TO CHECK OK
  Mmin=M-magbin/2.0;
  Mmax=M+magbin/2.0;
  Mstep=(Mmax-Mmin)/(npoints);
  nmean=0.0;
  for(i=0;i<=npoints;i++){
     //cout<<Mmin<<"\t"<<Mmax<<"\t"<<Mmin+i*Mstep<<endl;
     nmean+=Gal.luminosityFunctionM(Mmin+i*Mstep,alpha0,phi0,M0)*Mstep;
  }
  //cout<<M<<"\t"<<nmean<<endl;
  return nmean;
}

//window function for pencil beam survey
//perform the integration by gridding and summation in crude fashion
double FisherLF::getWindow(CosmParam *c)
{
   int i,j,l;
   double x,xmin, xmax, xstep;
   double y,ymin, ymax, ystep;
   double z,zmin, zmax, zstep;
   double wx,wy,wz;
   double k;
   double dx, dy, dz;
   double dvol, pk;
   double window(0.0),wijk;
   int npoints(40);

   Cosmology cosm(c->omm,c->oml,c->omb,c->h,c->sigma8,c->nscal,c->omnu);

   //cout<<"window"<<endl;
   //volume dimensions
   dx=coAngDiamDistance(exp_z,c)*anglex;
   dy=coAngDiamDistance(exp_z,c)*angley;
   dz=coordDistanceNum(max_z,c)-coordDistanceNum(min_z,c);
   cout<<"Dim: "<<dx<<"\t"<<dy<<"\t"<<dz<<endl;

   //grid in log space
   xmin=log(1.0e-3);
   ymin=log(1.0e-3);
   zmin=log(1.0e-3);
   xmax=log(1.0e1);   
   ymax=log(1.0e1);
   zmax=log(1.0e1);

   xstep=(xmax-xmin)/(double)(npoints-1);
   ystep=(ymax-ymin)/(double)(npoints-1);
   zstep=(zmax-zmin)/(double)(npoints-1);
   dvol=xstep*ystep*zstep;
   
   for(i=0;i<npoints;i++){
      for(j=0;j<npoints;j++){
         for(l=0;l<npoints;l++){
            x=exp(xmin+i*xstep);
            y=exp(ymin+j*ystep);
            z=exp(zmin+l*zstep);
            wx=2.0*sin(x*dx/2.0)/(x*dx);
            wy=2.0*sin(y*dy/2.0)/(y*dy);
            wz=2.0*sin(z*dz/2.0)/(z*dz);
            k=sqrt(x*x+y*y+z*z);
            pk=Plinear(k,z,&cosm);
            wijk=pow(wx*wy*wz,2.0);
            wijk*=x*y*k;  //log steps
            wijk*=pk;
            window+=wijk;
         }
      }
   }
   window/=pow(2.0*PI,3.0);
   window*=dvol;
   
   return window;
}
//////////////////////////////////////////////////////////////////////
// Tinker bias model
//////////////////////////////////////////////////////////////////////
//fit to bias from Tinker 2008
double FisherLF::biasTinker(double mass, double z, CosmParam *cp)
{
   double nu, bias,dummy;
   double deltasc,sig;
   double A(1.0), a(0.1325);
   double B(0.183), b(1.5);
   double C(0.265), c(2.4);
   Cosmology cosm(cp->omm,cp->oml,cp->omb,cp->h,cp->sigma8,cp->nscal,cp->omnu);
   deltasc = cosm.delCrit0(z)/cosm.growthFac(z);
   sig = sigm(mass,dummy,&cosm,0);
   nu=deltasc/sig;

   bias=1.0;
   bias-= A*pow(nu,a)/(pow(nu,a)+pow(deltasc,b));
   bias+=B*pow(nu,b);
   bias+=C*pow(nu,c);

   return bias;
}


//////////////////////////////////////////////////////////////////////
// Abundance matching
//////////////////////////////////////////////////////////////////////
//
// Calculate galaxy bias using abundance matching:
// set n(>mass)=n(>Lum)=n(<mag)
// calculating the LHS from the halo mass function and 
// the RHS from the luminosity function
//

//calculate bias of galaxies in magnitude bin M +- dM at redshift z
double FisherLF::biasM(double M, double dM, double z, CosmParam *c)
{
   double mp, mm;
   double bp, bm;
   double np, nm;

   double bias;
   Cosmology cosm(c->omm,c->oml,c->omb,c->h,c->sigma8,c->nscal,c->omnu);

   mp=getHaloMassFromMag(M+dM/2.0,z,c);
   mm=getHaloMassFromMag(M-dM/2.0,z,c);
   
   //cout<<"M: "<<mp<<"\t"<<mm<<endl;

   //in mean bias and number calc mass sets lower limit
   bp=meanHaloBiasNW(z,&cosm,mp);
   bm=meanHaloBiasNW(z,&cosm,mm);
   //cout<<"Bpm: "<<bp<<"\t"<<bm<<endl;

   //np=cosm.nCollObject(z,mp);
   //nm=cosm.nCollObject(z,mm);
   np=meanHaloNumber(z,&cosm,mp);
   nm=meanHaloNumber(z,&cosm,mm);


   bias=(bm-bp)/(nm-np);

   return bias;
}

//get mass corresponding to mag from abundance matching
double FisherLF::getHaloMassFromMag(double M, double z, CosmParam *c)
{
   double mass;
   Galaxies Gal;
   double nM;
   double lphi0,phi0, M0, alpha0;
   double m1,m2;
   double tol(1.0e-4);

  //get parameters from tags specified when reioniztion model is set
  lphi0=getParamFromTag("_lphi0",c);
  M0=getParamFromTag("_M0",c);
  alpha0=getParamFromTag("_alpha0",c);
  phi0=exp(lphi0*log(10.0));

  //cout<<M0<<"\t"<<alpha0<<"\t"<<phi0<<endl;

  Cosmology cosm(c->omm,c->oml,c->omb,c->h,c->sigma8,c->nscal,c->omnu);
  setDummyGetMass(z,&cosm,1);
  Gal.setCosmology(&cosm);

   nM=Gal.numberDensityLF_M(M,alpha0,phi0,M0);
   m1=coolMass(&cosm,z);
   m2=1.0e5*m1;
   while(dummyGetMass(m1)>nM){
      m1*=10.0;
      //cout<<m1<<"\t"<<m2<<"\t"<<dummyGetMass(m1)<<"\t"<<nM<<endl;
   }
   m1/=10.0;


   //cout<<nM<<"\t"<<dummyGetMass(m1)<<"\t"<<dummyGetMass(m2)<<endl;
   mass=zriddrConst(dummyGetMass,nM,m1,m2,tol);
   //cout<<mass<<endl;
   return mass;
}

double dummyGetMass(double mass)
{
   return setDummyGetMass(mass,NULL,0);
}

double setDummyGetMass(double mass, Cosmology *c1, int iflag)
{
   static Cosmology *c;
   static double z;
   double ntot;
   double ntot2;

   if(iflag==1){
      z=mass;
      c=c1;
      return 0.0;
   }

   //ntot=c->nCollObject(z,mass);  //only set up for PS mass fcn -slow but accurate
   ntot=meanHaloNumber(z,c,mass);
   //cout<<ntot<<"\t"<<ntot2<<endl;

   return ntot;
}
