
 /* 
  * MAke use of Galaxies class to plot some simple Luminosity functions
  * and observational constraints
  *
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
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
// Initialise some basic quantities
///////////////////////////////////////////////////////////////
    
     double phi0, M0, alpha0;
    stringstream ss;
    FisherLF lf;
  string dirbase;

  Galaxies Gal;

  
  dirbase="./lfdata";
  CosmParam fiducial;
  string namec="_fiducial";

  lf.setDirbase(dirbase);

  file="./project.ini";
  fiducial=lf.setCosmParamFromFile(file,namec);
  Cosmology cosm(fiducial.omm,fiducial.oml,fiducial.omb,fiducial.h,fiducial.sigma8,fiducial.nscal,fiducial.omnu);
  Gal.initGalaxies(&cosm,NULL);

  ///////////////////////////////////////////////////////////////////////
// Use Bowyens fits to LF parameters to output the high-z galaxy LF
/////////////////////////////////////////////////////////////////////////
     vector<double> LFparam;
     double z,zmin(6.0), dz(1.0);
     double mab,L0;
     double Lmin,Lmax;  //limits in luminosity
     double m, mstep, mmin, mmax; //limits in apparent magnitude
     int nm(35);
     double nL, L, M, nM;
     double LM;
     double vol;

     mmin=20.0;
     mmax=34.0;
     mstep=(mmax-mmin)/(double)(nm);

     for(i=0;i<=10;i++){
        z=zmin+dz*i;
        LFparam=Gal.getLFParam(z);
        M0=LFparam[0];
        phi0=LFparam[1];
        alpha0=LFparam[2];

        //convert to luminosity for LF
        L0=Gal.lumFromAbsMag(M0,z);

        ss.clear();
        ss<<dirbase<<"/lf_z"<<z<<".dat";
        ss>>file;
        fout.open(file.c_str());
        cout<<file<<endl;
        cout<<"# "<<z<<"\t"<<M0<<"\t"<<phi0<<"\t"<<alpha0<<"\t"<<L0<<endl;
        fout<<"# "<<z<<"\t"<<M0<<"\t"<<phi0<<"\t"<<alpha0<<"\t"<<L0<<endl;
        for(j=0;j<=nm;j++){

           m=mmin+j*mstep;
           L=Gal.lumFromMag(m,z);
           nL=Gal.luminosityFunction(L,alpha0,phi0,L0);

           mab=Gal.magFromL(L,z); //relative magnitude
           M=Gal.absMagFromL(L,z); //absolute magnitude
           nM=Gal.luminosityFunctionM(M,alpha0,phi0,M0);
           LM=Gal.lumFromAbsMag(M,z);
           //conversion to unit angle
           vol=cosm.volumeComoving(z)*pow(PI/180.0,2.0);
           cout<<"check L: "<<L<<"\t"<<LM<<"\t"<<nM<<"\t"<<nL<<endl;
           //fout<<L<<"\t"<<mab<<"\t"<<M<<"\t"<<log10(nL*vol)<<"\t"<<log10(nM*vol)<<endl;
           fout<<L<<"\t"<<mab<<"\t"<<M<<"\t"<<nL*vol<<"\t"<<nM*vol<<endl;
           cout<<"check: "<<m<<"\t"<<mab<<"\t"<<log10(nL*vol)<<"\t"<<log10(nM*vol)<<endl;

        }//end L loop
        fout.close();

     }  

///////////////////////////////////////////////////////////////////////
// Output regions constrained by different high-z galaxy surveys
/////////////////////////////////////////////////////////////////////////
/*
  double survey_area;
  double survey_depth;
  double survey_limit;

  vector<double> areas, depths;
  vector<string> names;
  string survey_name;

  //double m, mstep, mmin, mmax; //limits in apparent magnitude
  //int nm(35);
  
  nm=150;
  mmin=16.0;
  mmax=34.0;
  mstep=(mmax-mmin)/(double)(nm);

  z=6.0; //used for M to m conversions on GOODS and UDF depth - Should be fixed

  //specify survey area in deg^2 and survey depth in apparent magnitude
  //GOODS
  survey_name="goods";
  survey_area=2.0*(10.0/60.0)*(16.0/60.0);
  survey_depth=Gal.magFromL(Gal.lumFromAbsMag(-20.0,z),z);
  names.push_back(survey_name);
  areas.push_back(survey_area);
  depths.push_back(survey_depth);
  //UDF
  survey_name="udf";
  survey_area=11.0*(1.0/60.0)*(1.0/60.0);
  survey_depth=Gal.magFromL(Gal.lumFromAbsMag(-18.0,z),z);
  names.push_back(survey_name);
  areas.push_back(survey_area);
  depths.push_back(survey_depth);
  //JWST
  survey_name="jwst";
  survey_area=2.0*(10.0/60.0)*(16.0/60.0);
  survey_depth=31.0;
  names.push_back(survey_name);
  areas.push_back(survey_area);
  depths.push_back(survey_depth);
  //Euclid wide
  survey_name="euclid_wide";
  survey_area=20000.0;
  survey_depth=24.5;
  names.push_back(survey_name);
  areas.push_back(survey_area);
  depths.push_back(survey_depth);
  //Euclid deep
  survey_name="euclid_deep";
  survey_area=40.0;
  survey_depth=26.5;
  names.push_back(survey_name);
  areas.push_back(survey_area);
  depths.push_back(survey_depth);

  for(i=0;i<names.size();i++){
     ss.clear();
     ss<<dirbase<<"/"<<names[i]<<"_region.dat"<<endl;
     ss>>file;
     cout<<file<<endl;
     fout.open(file.c_str());
     survey_area=areas[i];
     survey_depth=depths[i];
     fout<<"# "<<survey_area<<"\t"<<survey_depth<<endl;
     for(j=0;j<=nm;j++){

        m=mmin+j*mstep;
        if(m>survey_depth){
           survey_limit=1.0e30;
        }else{
           survey_limit=1.0/survey_area;
        }
        fout<<m<<"\t"<<survey_limit<<endl;

     }
     fout.close();
  }
*/
///////////////////////////////////////////////////////////////////////
// Output Quasar number densities
/////////////////////////////////////////////////////////////////////////
   double Lsol(3.9e33);
        double lphi0, lL0, gamma1, gamma2; //QSO LF parameters
     //follow same proceedure as for galaxies
     mmin=20.0;
     mmax=26.0;
     mstep=(mmax-mmin)/(double)(nm);

     for(i=0;i<=10;i++){
        z=zmin+dz*i;
        LFparam=Gal.getLFParamQSO(z);

        lphi0=LFparam[0];
        lL0=LFparam[1];
        L0=pow(10.0,lL0)*Lsol;
        phi0=pow(10.0,lphi0);
        gamma1=LFparam[2];
        gamma2=LFparam[3];

        ss.clear();
        ss<<dirbase<<"/lf_qso_z"<<z<<".dat";
        ss>>file;
        fout.open(file.c_str());
        cout<<file<<endl;
        cout<<"# "<<z<<"\t"<<L0<<"\t"<<phi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<lL0<<"\t"<<lphi0<<endl;
        fout<<"# "<<z<<"\t"<<L0<<"\t"<<phi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<lL0<<"\t"<<lphi0<<endl;
        for(j=0;j<=nm;j++){

           m=mmin+j*mstep;
           L=Gal.lumFromMag(m,z);
           nL=Gal.luminosityFunctionQSO_Hop(L,phi0,L0,gamma1,gamma2);

           mab=Gal.magFromL(L,z); //relative magnitude
           M=Gal.absMagFromL(L,z); //absolute magnitude

           //need to convert dNdL to dNdM
           //nM=Gal.luminosityFunctionM(M,alpha0,phi0,M0);
           nM=nL*(1.08574/L);

           //cout<<nM<<"\t"<<nL<<endl;

           LM=Gal.lumFromAbsMag(M,z);
           //conversion to unit angle
           vol=cosm.volumeComoving(z)*pow(PI/180.0,2.0);
           //fout<<L<<"\t"<<mab<<"\t"<<M<<"\t"<<log10(nL*vol)<<"\t"<<log10(nM*vol)<<endl;
           fout<<L<<"\t"<<mab<<"\t"<<M<<"\t"<<nL*vol<<"\t"<<nM*vol<<endl;
           //cout<<"check: "<<m<<"\t"<<mab<<"\t"<<log10(nL*vol)<<"\t"<<log10(nM*vol)<<endl;

        }//end L loop
        fout.close();

     } 



  ///////////////////////////////////////////////////////////////////////
  return 0;
}




