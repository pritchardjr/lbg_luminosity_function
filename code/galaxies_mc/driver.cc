/* driver.cc
 * 
 * Program for computing the fisher matrix for a galaxy survey.
 * Call by
 *   ./driver.x -z6 > (output filename)
 */

#include <iostream>
#include <fstream.h>
#include "dnumrecipes.h"
#include "math.h"
#include "dcosmology.h"
#include "fisherCMB.h"
#include "fisherGAL.h"
#include "fisher21CM.h"
#include "fisher.h"
#include "bubble.h"
#include "haloDensity.h"

int bubble_flag(0);
int bubble_use(0);

using namespace std;

const double INT_ACCURACY(1.0e-6);

int main(int argc, char *argv[])
{
  double z;
  int i;
  int j;
  int reset_flag(0);
  char *file;
  ofstream fout;
  /////////////////////////////////////////////////////////////
  // Handle command line input
  ////////////////////////////////////////////////////////////

  while ((argc>1) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'r':
      reset_flag = atoi(&argv[1][2]);
      break;
    case 'e':
      //reset_flag = atoi(&argv[1][2]);
      reset_flag=2;
      break;
    case 'b':
      //reset_flag = atoi(&argv[1][2]);
      bubble_flag=1;
      bubble_use=1;
      cout<<"Using bubbles"<<endl;
      break;
    case 'w':
      //use bubble in ps but not fisher analysis
      bubble_use=1;
      cout<<"Using bubbles in power spectrum"<<endl;
      break;
    default:
      cerr << "Bad option " << argv[1] << "\n";
    }
    --argc;
    ++argv;
  }

  /////////////////////////////////////////////////////////////////
  // Specify fiducial cosmology
  /////////////////////////////////////////////////////////////////
  FisherCMB WMAPF;
  FisherCMB PLANCK;
  string name;

  CosmParam fiducial;
  double h,omm,omb,oml,omnu,tau,ascal2,nscal,alpha,ts,w0,w1,sigma8;
  double Mnu;
  string namec="_fiducial";

  h=0.7;
  omm=0.3;
  omb=0.046;
  oml=0.7;
  Mnu=0.3;
  omnu=WMAPF.omnuhhFromMnu(Mnu)/h/h;
  tau=0.1;
  nscal=0.95;
  alpha=0.0;
  ascal2=29.5;
  ascal2*=0.9;
  ts=0.0;
  w0=-1.0;
  w1=0.0;

/*
//lesgourgues parameters
  h=0.69;
  omm=0.3;
  omb=0.048;
  oml=0.7;
  Mnu=0.3;
//  Mnu=0.1;
  omnu=WMAPF.omnuhhFromMnu(Mnu)/h/h;
  tau=0.11;
  nscal=0.96;
  alpha=0.0;
  ascal2=29.5;
  ascal2*=0.9;
  ts=0.0;
  w0=-1.0;
  w1=0.0;
*/

if(reset_flag==1){
  fiducial=WMAPF.setCosmParam(h,omm,omb,oml,omnu,tau,ascal2,nscal,alpha,ts,w0,w1,-1.0,namec);
}else{
 fiducial=WMAPF.setCosmParam(h,omm,omb,oml,omnu,tau,ascal2,nscal,alpha,ts,w0,w1,0.7693057,namec);
}

  Cosmology cosm(omm,oml,omb,h,fiducial.sigma8,nscal,omnu);

// specify parameters
  int *choiceVec;
  choiceVec=ivector(1,22);
  for(i=1;i<=22;i++) choiceVec[i]=1;
  choiceVec[6]=1;   //cosmological constant
  choiceVec[7]=1;   // tau
  choiceVec[8]=0;   //tensors
  choiceVec[9]=0;   //flat Universe
  choiceVec[10]=1;   // helium
  choiceVec[11]=1;  //running
  choiceVec[12]=1;  //Mnu
  choiceVec[13]=0;  // bias
  choiceVec[14]=0;  //numass1
  choiceVec[15]=0;  //numass2
  choiceVec[16]=0;  //numass3
  choiceVec[17]=0;  //no lda
  choiceVec[18]=0;  //no lh
  choiceVec[19]=0;  //no lg   
  choiceVec[20]=0;  //epsilon 
  choiceVec[21]=0;  //eta
  choiceVec[22]=0;  //xi 

  WMAPF.specifyParams(choiceVec,fiducial);

  cout<<"initialising experiments"<<endl;
  //initialise different cosmologies for derivatives
  if(reset_flag==1) WMAPF.setCosmFiles();

  /////////////////////////////////////////////////////////////////
  // CMB Experiment Parameters
  ////////////////////////////////////////////////////////////////
  choiceVec[6]=1;   //cosmological constant
  choiceVec[11]=0;  //running
  choiceVec[12]=1;  //Mnu
  choiceVec[14]=0;  //numass1
  choiceVec[15]=0;  //numass2
  choiceVec[16]=0;  //numass3

  int lmax(2500);

  cout<<"calculating Fisher matrix"<<endl;
//  WMAPF.fisherMatrixCMB();

  PLANCK.specifyParams(choiceVec,fiducial);
  name="./data/planck";
  PLANCK.initExperiment(name,0.65,3,lmax,143,7.1,6.0,11.4);
  PLANCK.addChannelExperiment(217,5.0,13.1,26.7);
  PLANCK.addChannelExperiment(100,9.5,6.8,10.9);

  PLANCK.fisherMatrixCMB();

  PLANCK.printErrors();

  PLANCK.latexTable();

  FisherCMB CosmicVar;

  CosmicVar.specifyParams(choiceVec,fiducial);
  name="./data/cosmicvar";

  CosmicVar.initExperiment(name,1.0,3,lmax,30,1.0e-10,1.0e-10,1.0e-10);

  CosmicVar.fisherMatrixCMB();

  CosmicVar.printErrors();

  CosmicVar.latexTable();

  /////////////////////////////////////////////////////////////////
  // Galaxy Experiment Parameters
  ////////////////////////////////////////////////////////////////
  choiceVec[13]=1;  //bias
  choiceVec[10]=0;   // no helium
  choiceVec[7]=0;   // no tau
  choiceVec[11]=0;  //no running
  choiceVec[5]=1;   //ascal2

  FisherGAL SDSS;
  double volume, density, kmin, kmax, sigma8gal,zuse, sigmaz;

  sigmaz=0.0;

  name="./data/sdss";
  volume=1.0e9*pow(fiducial.h,-3.0);
  density=1.0e-4*pow(fiducial.h,3.0);
  kmin=0.0003;
  kmax=0.1*fiducial.h;
  sigma8gal=1.8;
  zuse=0.3;
  SDSS.initSurvey(name, zuse, volume, density, kmin, kmax, sigma8gal, sigmaz);

	SDSS.specifyParams(choiceVec,fiducial);

	 SDSS.addParam("_bias");
	 SDSS.removeParam("_tau");
	 SDSS.removeParam("_ts");
	SDSS.fisherMatrixGAL();

  SDSS.printErrors();
  SDSS.latexTable();
  SDSS.outPowerErrors();

  Fisher SDSS_PLANCK;
  name="./data/sdss_planck";
  SDSS_PLANCK.setName(name);

  SDSS_PLANCK.combineFisherMatrices(SDSS.returnFisherData(), PLANCK.returnFisherData(), &fiducial);
  SDSS_PLANCK.latexTable();


  FisherGAL LSST;

  name="./data/lsst";
  volume=1.8e9*pow(fiducial.h,-3.0);
  density=5.0e-4*pow(fiducial.h,3.0);
  kmin=0.0003;
  kmax=0.2*fiducial.h;
  sigma8gal=1.0;
  zuse=1.0;
  LSST.initSurvey(name, zuse, volume, density, kmin, kmax, sigma8gal, sigmaz);

	LSST.specifyParams(choiceVec,fiducial);

	 LSST.addParam("_bias");
	 LSST.removeParam("_tau");
	 LSST.removeParam("_ts");
	LSST.fisherMatrixGAL();

  LSST.printErrors();
  LSST.latexTable();
  LSST.outPowerErrors();


  Fisher LSST_PLANCK;
  name="./data/lsst_planck";
  LSST_PLANCK.setName(name);

  LSST_PLANCK.combineFisherMatrices(LSST.returnFisherData(), PLANCK.returnFisherData(), &fiducial);
  LSST_PLANCK.latexTable();


  FisherGAL LYB;

  name="./data/lyb";
  volume=0.5e9*pow(fiducial.h,-3.0);
  density=1.0e-3*pow(fiducial.h,3.0);
  kmin=0.0003;
  kmax=0.5*fiducial.h;
  sigma8gal=1.0;
  zuse=3.0;
  LYB.initSurvey(name, zuse, volume, density, kmin, kmax, sigma8gal, sigmaz);

	LYB.specifyParams(choiceVec,fiducial);

	 LYB.addParam("_bias");
	 LYB.removeParam("_tau");
	 LYB.removeParam("_ts");
	LYB.fisherMatrixGAL();

  LYB.printErrors();
  LYB.latexTable();
  LYB.outPowerErrors();


  Fisher LYB_PLANCK;
  name="./data/lyb_planck";
  LYB_PLANCK.setName(name);

  LYB_PLANCK.combineFisherMatrices(LYB.returnFisherData(), PLANCK.returnFisherData(), &fiducial);
  LYB_PLANCK.latexTable();

  FisherGAL SUPERGAL;

  name="./data/supergal";
  volume=40.0e9*pow(fiducial.h,-3.0);
  density=5.0e-4*pow(fiducial.h,3.0);
  kmin=0.0003;
  kmax=0.2*fiducial.h;
  sigma8gal=1.0;
  zuse=1.0;
  SUPERGAL.initSurvey(name, zuse, volume, density, kmin, kmax, sigma8gal, sigmaz);

	SUPERGAL.specifyParams(choiceVec,fiducial);

	 SUPERGAL.addParam("_bias");
	 SUPERGAL.removeParam("_tau");
	 SUPERGAL.removeParam("_ts");
	SUPERGAL.fisherMatrixGAL();

  SUPERGAL.printErrors();
  SUPERGAL.latexTable();
  SUPERGAL.outPowerErrors();

  Fisher SUPERGAL_PLANCK;
  name="./data/supergal_planck";
  SUPERGAL_PLANCK.setName(name);

  SUPERGAL_PLANCK.combineFisherMatrices(SUPERGAL.returnFisherData(), PLANCK.returnFisherData(), &fiducial);
  SUPERGAL_PLANCK.latexTable();


  FisherGAL HIGHGAL;

  name="./data/highgal";
  volume=20.0e9*pow(fiducial.h,-3.0);
  density=1.0e-3*pow(fiducial.h,3.0);
  kmin=0.0003;
  kmax=0.35*fiducial.h;
  sigma8gal=1.0;
  zuse=3.0;
  HIGHGAL.initSurvey(name, zuse, volume, density, kmin, kmax, sigma8gal, sigmaz);

	HIGHGAL.specifyParams(choiceVec,fiducial);

	 HIGHGAL.addParam("_bias");
	 HIGHGAL.removeParam("_tau");
	 HIGHGAL.removeParam("_ts");
	HIGHGAL.fisherMatrixGAL();

  HIGHGAL.printErrors();
  HIGHGAL.latexTable();
  HIGHGAL.outPowerErrors();

  Fisher HIGHGAL_PLANCK;
  name="./data/highgal_planck";
  HIGHGAL_PLANCK.setName(name);

  HIGHGAL_PLANCK.combineFisherMatrices(HIGHGAL.returnFisherData(), PLANCK.returnFisherData(), &fiducial);
  HIGHGAL_PLANCK.latexTable();

//////////////////////////////////////////////////////////////////
// 21 cm experiments
//////////////////////////////////////////////////////////////////

 string files;


 Fisher21CM MWA;
 Observation ObsMWA(&cosm);

 files="mwa_param.ini";
 ObsMWA.initObservationFromFile(files);

//  double kmin,kmax;
  kmin=0.0003;
  kmax=2.0; // *fiducial.h;

 name="./data/mwa";
 MWA.initObservation(name,8.0,kmin,kmax,&ObsMWA,0);

	MWA.specifyParams(choiceVec,fiducial);

	 MWA.addParam("_bias");
	 MWA.removeParam("_tau");
	 MWA.removeParam("_ts");
	MWA.fisherMatrix21CM();

	MWA.latexTable();
	MWA.outPowerErrors();



 Fisher21CM SKA;
 Observation ObsSKA(&cosm);

 files="ska_param.ini";
 ObsSKA.initObservationFromFile(files);

  kmin=0.0003;
  kmax=2.0; // *fiducial.h;

 name="./data/ska";
 SKA.initObservation(name,8.0,kmin,kmax,&ObsSKA,0);

	SKA.specifyParams(choiceVec,fiducial);

	 SKA.addParam("_bias");
	 SKA.removeParam("_tau");
	 SKA.removeParam("_ts");
	SKA.fisherMatrix21CM();

	SKA.latexTable();
	SKA.outPowerErrors();

 Fisher21CM LOFAR;
 Observation ObsLOFAR(&cosm);

 files="lofar_param.ini";
 ObsLOFAR.initObservationFromFile(files);

  kmin=0.0003;
  kmax=2.0; // *fiducial.h;

 name="./data/lofar";
 LOFAR.initObservation(name,8.0,kmin,kmax,&ObsLOFAR,0);

	LOFAR.specifyParams(choiceVec,fiducial);

	 LOFAR.addParam("_bias");
	 LOFAR.removeParam("_tau");
	 LOFAR.removeParam("_ts");
	LOFAR.fisherMatrix21CM();

	LOFAR.latexTable();
	LOFAR.outPowerErrors();


 Fisher21CM LA;
 Observation ObsLA(&cosm);

 files="la_param.ini";
 ObsLA.initObservationFromFile(files);

  kmin=0.0003;
  kmax=2.0; // *fiducial.h;

 name="./data/la";
 LA.initObservation(name,8.0,kmin,kmax,&ObsLA,0);

	LA.specifyParams(choiceVec,fiducial);

	 LA.addParam("_bias");
	 LA.removeParam("_tau");
	 LA.removeParam("_ts");
	LA.fisherMatrix21CM();

	LA.latexTable();
	LA.outPowerErrors();

  Fisher LA_PLANCK;
 name="./data/la_planck";
  LA_PLANCK.setName(name);

  LA_PLANCK.combineFisherMatrices(LA.returnFisherData(), PLANCK.returnFisherData(), &fiducial);

  LA_PLANCK.printErrors();

  LA_PLANCK.latexTable();

///////////////////////////////////////////////////////////////////
// Combine fisher matrices
///////////////////////////////////////////////////////////////////

  Fisher MWA_PLANCK;

 name="./data/mwa_planck";
  MWA_PLANCK.setName(name);

  MWA_PLANCK.combineFisherMatrices(MWA.returnFisherData(), PLANCK.returnFisherData(),&fiducial);

  MWA_PLANCK.printErrors();
  MWA_PLANCK.latexTable();


  Fisher SKA_PLANCK;
 name="./data/ska_planck";
  SKA_PLANCK.setName(name);

  SKA_PLANCK.combineFisherMatrices(SKA.returnFisherData(), PLANCK.returnFisherData(), &fiducial);

  SKA_PLANCK.printErrors();
  SKA_PLANCK.latexTable();
///////////////////////////////////////////////////////////////////
// SKA redshift slices
///////////////////////////////////////////////////////////////////
//string files;
//double kmin, kmax;

  kmin=0.0003;
  kmax=2.0; // *fiducial.h;


// Observation ObsSKA(&cosm);

 files="ska_param.ini";
 //ObsSKA.initObservationFromFile(files);

 Fisher21CM SKA1, SKA2, SKA3;

  //specify parameters
//  for(i=1;i<=13;i++) choiceVec[i]=1;
//  choiceVec[6]=0;   //cosmological constant
//  choiceVec[8]=0;   //no tensors
//  choiceVec[9]=0;   //flat Universe
//  choiceVec[7]=0;   // no tau

 name="./data/ska1";
 SKA1.initObservation(name,8.0,kmin,kmax,&ObsSKA,0);
 SKA1.specifyParams(choiceVec,fiducial);
 SKA1.fisherMatrix21CM();
 SKA1.outPowerErrors();

 name="./data/ska2";
 SKA2.initObservation(name,7.5,kmin,kmax,&ObsSKA,0);
 SKA2.specifyParams(choiceVec,fiducial);
 SKA2.fisherMatrix21CM();
 SKA2.outPowerErrors();

 name="./data/ska3";
 SKA3.initObservation(name,7.0,kmin,kmax,&ObsSKA,0);
 SKA3.specifyParams(choiceVec,fiducial);
 SKA3.fisherMatrix21CM();
 SKA3.outPowerErrors();

  Fisher SKA_FULL;
  name="./data/ska_full";
  SKA_FULL.setName(name);

  SKA_FULL.combineFisherMatrices(SKA1.returnFisherData(), SKA2.returnFisherData(), &fiducial);
  SKA_FULL.combineFisherMatrices(SKA_FULL.returnFisherData(), SKA3.returnFisherData(), &fiducial);

  SKA_FULL.printErrors();

  Fisher SKA_FULL_PLANCK;
 name="./data/ska_full_planck";
  SKA_FULL_PLANCK.setName(name);

  SKA_FULL_PLANCK.combineFisherMatrices(SKA_FULL.returnFisherData(), PLANCK.returnFisherData(), &fiducial);

  SKA_FULL_PLANCK.printErrors();

  SKA_FULL.latexTable();
  SKA_FULL_PLANCK.latexTable();

///////////////////////////////////////////////////////////////////
// FFTT redshift slices
///////////////////////////////////////////////////////////////////
//string files;
//double kmin, kmax;

//  kmin=0.0003;
//  kmax=2.0; // *fiducial.h;

 Observation ObsFFTT(&cosm);

 files="fftt_param.ini";
 ObsFFTT.initObservationFromFile(files);

 Fisher21CM FFTT1, FFTT2, FFTT3;

 name="./data/fftt1";
 FFTT1.initObservation(name,8.0,kmin,kmax,&ObsFFTT,0);
 FFTT1.specifyParams(choiceVec,fiducial);
 FFTT1.fisherMatrix21CM();
 FFTT1.outPowerErrors();

 FFTT1.latexTable();

  Fisher FFTT1_PLANCK;
 name="./data/fftt1_planck";
  FFTT1_PLANCK.setName(name);

  FFTT1_PLANCK.combineFisherMatrices(FFTT1.returnFisherData(), PLANCK.returnFisherData(), &fiducial);

  FFTT1_PLANCK.printErrors();
  FFTT1_PLANCK.latexTable();	


 name="./data/fftt2";
 FFTT2.initObservation(name,7.5,kmin,kmax,&ObsFFTT,0);
 FFTT2.specifyParams(choiceVec,fiducial);
 FFTT2.fisherMatrix21CM();
 FFTT2.outPowerErrors();

 name="./data/fftt3";
 FFTT3.initObservation(name,7.0,kmin,kmax,&ObsFFTT,0);
 FFTT3.specifyParams(choiceVec,fiducial);
 FFTT3.fisherMatrix21CM();
 FFTT3.outPowerErrors();

  Fisher FFTT_FULL;
  name="./data/fftt_full";
  FFTT_FULL.setName(name);

  FFTT_FULL.combineFisherMatrices(FFTT1.returnFisherData(), FFTT2.returnFisherData(), &fiducial);
  FFTT_FULL.combineFisherMatrices(FFTT_FULL.returnFisherData(), FFTT3.returnFisherData(), &fiducial);

  FFTT_FULL.printErrors();

  Fisher FFTT_FULL_PLANCK;
 name="./data/fftt_full_planck";
  FFTT_FULL_PLANCK.setName(name);

  FFTT_FULL_PLANCK.combineFisherMatrices(FFTT_FULL.returnFisherData(), PLANCK.returnFisherData(), &fiducial);

  FFTT_FULL_PLANCK.printErrors();

  FFTT_FULL.latexTable();
  FFTT_FULL_PLANCK.latexTable();

  Fisher FFTT_FULL_CosmicVar;
 name="./data/fftt_full_cosmicvar";
  FFTT_FULL_CosmicVar.setName(name);

  FFTT_FULL_CosmicVar.combineFisherMatrices(FFTT_FULL.returnFisherData(), CosmicVar.returnFisherData(), &fiducial);

  FFTT_FULL_CosmicVar.printErrors();

  FFTT_FULL.latexTable();
  FFTT_FULL_CosmicVar.latexTable();

///////////////////////////////////////////////////////////////////
// FFTT_HIGH redshift slices
///////////////////////////////////////////////////////////////////
 /*
 string dirbase_inf;
 Observation ObsFFTT_HIGH(&cosm);

 dirbase_inf="./data";

 files="fftt_high_param.ini";
 ObsFFTT_HIGH.initObservationFromFile(files);

 Fisher21CM FFTT_HIGH1, FFTT_HIGH2, FFTT_HIGH3;

  FFTT_HIGH1.setDirbase(dirbase_inf);
  FFTT_HIGH2.setDirbase(dirbase_inf);
  FFTT_HIGH3.setDirbase(dirbase_inf);

  kmin=0.0003;
  kmax=40.0; // *fiducial.h;
  zuse=20.5;
 name=dirbase_inf+"/fftt_high1";
 FFTT_HIGH1.initObservation(name,zuse,kmin,kmax,&ObsFFTT_HIGH,0);
 FFTT_HIGH1.specifyParams(choiceVec,fiducial);
 FFTT_HIGH1.fisherMatrix21CM();
 FFTT_HIGH1.outPowerErrors();

 FFTT_HIGH1.latexTable();

  Fisher FFTT_HIGH1_PLANCK;
 name=dirbase_inf+"/fftt_high1_planck";
  FFTT_HIGH1_PLANCK.setName(name);

  FFTT_HIGH1_PLANCK.combineFisherMatrices(FFTT_HIGH1.returnFisherData(), PLANCK.returnFisherData(), &fiducial);

  FFTT_HIGH1_PLANCK.printErrors();
  FFTT_HIGH1_PLANCK.latexTable();	

  kmin=0.0003;
  kmax=40.0; // *fiducial.h;
  zuse=18.0;
 name=dirbase_inf+"/fftt_high2";
 FFTT_HIGH2.initObservation(name,zuse,kmin,kmax,&ObsFFTT_HIGH,0);
 FFTT_HIGH2.specifyParams(choiceVec,fiducial);
 FFTT_HIGH2.fisherMatrix21CM();
 FFTT_HIGH2.outPowerErrors();

  kmin=0.0003;
  kmax=40.0; // *fiducial.h;
  zuse=16.0;
 name=dirbase_inf+"/fftt_high3";
 FFTT_HIGH3.initObservation(name,zuse,kmin,kmax,&ObsFFTT_HIGH,0);
 FFTT_HIGH3.specifyParams(choiceVec,fiducial);
 FFTT_HIGH3.fisherMatrix21CM();
 FFTT_HIGH3.outPowerErrors();

  Fisher FFTT_HIGH_FULL;
  name=dirbase_inf+"/fftt_high_full";
  FFTT_HIGH_FULL.setName(name);

  FFTT_HIGH_FULL.combineFisherMatrices(FFTT_HIGH1.returnFisherData(), FFTT_HIGH2.returnFisherData(), &fiducial);
  FFTT_HIGH_FULL.combineFisherMatrices(FFTT_HIGH_FULL.returnFisherData(), FFTT_HIGH3.returnFisherData(), &fiducial);

  FFTT_HIGH_FULL.printErrors();

  Fisher FFTT_HIGH_FULL_PLANCK;
 name=dirbase_inf+"/fftt_high_full_planck";
  FFTT_HIGH_FULL_PLANCK.setName(name);

  FFTT_HIGH_FULL_PLANCK.combineFisherMatrices(FFTT_HIGH_FULL.returnFisherData(), PLANCK.returnFisherData(), &fiducial);

  FFTT_HIGH_FULL_PLANCK.printErrors();

  FFTT_HIGH_FULL.latexTable();
  FFTT_HIGH_FULL_PLANCK.latexTable();

*/
///////////////////////////////////////////////////////////////
// Combine GAL +21CM + CMB
//////////////////////////////////////////////////////////////

Fisher SDSS_SKA_PLANCK;
 name="./data/sdss_ska_planck";
  SDSS_SKA_PLANCK.setName(name);

SDSS_SKA_PLANCK.combineFisherMatrices(SKA_FULL.returnFisherData(), SDSS_PLANCK.returnFisherData(), &fiducial);
SDSS_SKA_PLANCK.latexTable();


Fisher SDSS_FFTT_PLANCK;
 name="./data/sdss_fftt_planck";
  SDSS_FFTT_PLANCK.setName(name);

SDSS_FFTT_PLANCK.combineFisherMatrices(FFTT_FULL.returnFisherData(), SDSS_PLANCK.returnFisherData(), &fiducial);
SDSS_FFTT_PLANCK.latexTable();


Fisher SUPERGAL_SKA_PLANCK;
 name="./data/supergal_ska_planck";
  SUPERGAL_SKA_PLANCK.setName(name);

SUPERGAL_SKA_PLANCK.combineFisherMatrices(SKA_FULL.returnFisherData(), SUPERGAL_PLANCK.returnFisherData(), &fiducial);
SUPERGAL_SKA_PLANCK.latexTable();


Fisher SUPERGAL_FFTT_PLANCK;
 name="./data/supergal_fftt_planck";
  SUPERGAL_FFTT_PLANCK.setName(name);

SUPERGAL_FFTT_PLANCK.combineFisherMatrices(FFTT_FULL.returnFisherData(), SUPERGAL_PLANCK.returnFisherData(), &fiducial);
SUPERGAL_FFTT_PLANCK.latexTable();

////////////////////////////////////////////////////////////////
  //clean up memory and end 
  return 0;
}
