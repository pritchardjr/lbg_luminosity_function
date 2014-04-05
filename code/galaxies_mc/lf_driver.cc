
 /* driver.cc
 * program for mimicing the results of BL04
 * Call by 
 *   driver.x -z6 > (output filename)
 *
 */

#include <iostream>
#include <fstream>
#include "math.h"
#include "astrophysics.h"
#include "dnumrecipes.h"
#include "dcosmology.h"
#include "twentyonecm.h"
#include "spline.h"
#include <sstream>
#include "fisherLF.h"
#include "galaxies.h"
#include <gsl/gsl_sf.h>

using namespace std;

/******************************************************************/


int main(int argc, char *argv[])
{
 
  // Cosmology parameters

  string name;
  double zin_arg(15.0);
  int i,j;
  string file;
  ofstream fout;
  ifstream fin;
  int xin(1);
  int lyaxray_in(0);
  int popflag_in(0);
  int file_flag(0);
  int reset_flag(0);
  int case_flag(0);

  // Handle arguments 
  while ((argc>1) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'z':
      zin_arg=atof(&argv[1][2]);
      break;
    case 'r':
      reset_flag=atoi(&argv[1][2]);
      break;
    case 'c':
      case_flag=atoi(&argv[1][2]);
      break;
    case 'l':
      lyaxray_in=atoi(&argv[1][2]);
      break;
    case 'p':
      popflag_in=atoi(&argv[1][2]);
      break;
    case 'f':
      file_flag=1;
      break;
    default:
      cerr << "Bad option" <<argv[1] <<"\n";
    }
    --argc;
    ++argv;
  }

  

///////////////////////////////////////////////////////////////
// Calculate Fisher matrix constraints on LF for a simple survey
///////////////////////////////////////////////////////////////
    
     double lphi0,phi0, M0, alpha0;
    stringstream ss;
    double volume, maglim;
    double anglex, angley;
    int nfield;
    int nbin;
    FisherLF lf;
  string dirbase_inf;
  double zmin, zmax, area;
  double m,z;
  Galaxies Gal;

  
  dirbase_inf="./data_lf";
  CosmParam fiducial;
  string namec="_fiducial";

  lf.setDirbase(dirbase_inf);

  file="./project.ini";
  fiducial=lf.setCosmParamFromFile(file,namec);
  Cosmology cosm(fiducial.omm,fiducial.oml,fiducial.omb,fiducial.h,fiducial.sigma8,fiducial.nscal,fiducial.omnu);
  Gal.initGalaxies(&cosm,NULL);

// specify parameters for reset
  int *choiceVec;
  choiceVec=ivector(1,24);
  for(i=1;i<=24;i++) choiceVec[i]=0;

  lf.specifyParams(choiceVec,fiducial);

  lphi0=log10(1.0e-3);
  M0=-20;
  alpha0=-1.74;

  lf.setLFParam(lphi0,M0,alpha0);

  if(case_flag==1){
//UDF

  zmin=6.5;
  zmax=7.5;
  z=(zmax+zmin)/2.0;
  area=(4.7/60.0/60.0);
  anglex=2.2;
  angley=2.2;
  volume=lf.volumeFromArea(zmin,zmax,area);
  cout<<volume<<endl;
  m=29.14;
  maglim=Gal.absMagFromM(m,z);
  nfield=1;
  nbin=20;

    name=dirbase_inf+"/udf";
    lf.initExperiment(name,zmin,zmax,anglex,angley,maglim,nfield,nbin);
    lf.fisherMatrix();
    lf.latexTable();

  }else{
//GOODS

  zmin=6.5;
  zmax=7.5;
  z=(zmax+zmin)/2.0;
  area=(53.3/60.0/60.0);
  anglex=10.0;
  angley=5.33;
  volume=lf.volumeFromArea(zmin,zmax,area);
  cout<<volume<<endl;
  m=27.73;
  maglim=Gal.absMagFromM(m,z);
  nfield=1;
  nbin=20;

    name=dirbase_inf+"/goods";
    lf.initExperiment(name,zmin,zmax,anglex,angley,maglim,nfield,nbin);
    lf.fisherMatrix();
    lf.latexTable();
  }

  ///////////////////////////////////////////////////////////////////////
  return 0;
}




