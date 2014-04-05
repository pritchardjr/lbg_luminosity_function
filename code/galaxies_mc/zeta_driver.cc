
 /* driver.cc
 * program for mimicing the results of BL04
 * Call by 
 *   driver.x -z6 > (output filename)
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
#include "astrophysics.h"
#include "dnumrecipes.h"
#include "dcosmology.h"
#include "twentyonecm.h"
#include "reionization.h"
#include "haloDensity.h"
#include "spline.h"
#include "radio.h"
#include "observation.h"
#include "fisher.h"
#include "fisherGAL.h"
#include "neutrinos.h"
#include "lymanforest.h"
#include "ionization.h"
#include "spline2D.h"

using namespace std;

/******************************************************************/


int main(int argc, char *argv[])
{
 
  // Cosmology parameters
  double om0(0.3);
  double lam0(0.7);
  double omb(0.046);
  double h(0.7);
  double s8(0.9);
  double n(1.0);
  double omNu(0.0);

  double zin_arg(15.0);
  int tin(2);
  int i,j;
  char *file;
  ofstream fout;
  ifstream fin;
  int xin(1);
  int lyaxray_in(0);
  int popflag_in(0);

     int ndot_flag(0), cmb_flag(0);
     int param_flag(0);
     int lls_flag(1);
     int like_calc_flag(1);
     int nthread(0);
     int gamma_flag(0);
     int tocm_flag(0);

  // Handle arguments 
  while ((argc>1) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'z':
      zin_arg=atof(&argv[1][2]);
      break;
    case 't':
      nthread=atoi(&argv[1][2]);
      cout<<"thread no: "<<nthread<<endl;
      break;
    case 'x':
      xin=atoi(&argv[1][2]);
      break;
    case 'l':
      lls_flag=atoi(&argv[1][2]);
      break;
    case 'g':
      gamma_flag=atoi(&argv[1][2]);
      break;
    case 'p':
      param_flag=atoi(&argv[1][2]);
      break;
    case 'n':
      ndot_flag=atoi(&argv[1][2]);
      break;
    case 'c':
      cmb_flag=atoi(&argv[1][2]);
      break;
    case 'd':
      like_calc_flag=atoi(&argv[1][2]);
      break;
    case 'm':
      tocm_flag=atoi(&argv[1][2]);
      break;
    default:
      cerr << "Bad option" <<argv[1] <<"\n";
    }
    --argc;
    ++argv;
  }

  
  //  s8=0.77;
  //omb=0.044;
  //om0=0.26;
  //lam0=0.74;
  //h=0.72;
  //n=0.95;

  //HQB cosmology
  //  s8=0.9;
  //omb=0.044;
  //om0=0.27;
  //lam0=0.73;
  //h=0.71;
  //n=1.0;

  //NB07 cosmology
  //  s8=0.826;
  //omb=0.0478;
  //om0=0.299;
  //lam0=0.701;
  //h=0.687;
  //n=1.0;

  //Oh1999
  //  s8=0.87;
  //omb=0.04;
  //om0=0.35;
  //lam0=0.65;
  //h=0.65;
  //n=0.96;

  cout<<"Enter cosmology parameters"<<endl;
  cout<<"Enter s8:"<<endl;
  cin>>s8;
  cout<<"Enter omb:"<<endl;
  cin>>omb;
  cout<<"Enter om0:"<<endl;
  cin>>om0;
  cout<<"Enter lam0:"<<endl;
  cin>>lam0;
  cout<<"Enter h:"<<endl;
  cin>>h;
  cout<<"Enter n:"<<endl;
  cin>>n;

  cout<<s8<<"\t"<<omb<<"\t"<<om0<<"\t"<<lam0<<"\t"<<h<<"\t"<<n<<endl;
  
  Cosmology c(om0,lam0,omb,h,s8,n,omNu);
  Astrophysics a(&c,popflag_in,xin,lyaxray_in,1.0);
  TwentyOneCM tocm(&c,&a);
  Reionization rz(&c,&a);
  Lyman lyf;

  int lyaxrayflag, popflag, xrayflag, sourceflag;
  double fstar, fx, nion, nlya, fesc;
  int modelflag;

  modelflag=tin;
  if(modelflag==1){
    //Low tau parameters
    fstar=0.1;
    fesc=0.075;
    nion=4000.0;
    nlya=-1.0;
    fx=1.0;
    lyaxrayflag=0;
    //  lyaxrayflag=0;
    popflag=0; //popII
  }else if(modelflag==2){
    //Mid tau parameters
    fstar=0.2;
    fesc=0.2;
    nion=4000.0;
    nlya=-1.0;
    fx=1.0;
    lyaxrayflag=0;
    //  lyaxrayflag=0;
    popflag=0; //popIII
  }else if(modelflag==3){
    //High tau parameters
    fstar=0.3;
    fesc=0.1;
    nion=30000.0;
    nlya=3030.0;
    fx=1.0;
    lyaxrayflag=0;
    //  lyaxrayflag=0;
    popflag=1; //popIII
  }

  xrayflag=-1;
  sourceflag=-1;

  cout<<"Enter astrophysics parameters"<<endl;
  cout<<"Enter fstar:"<<endl;
  cin>>fstar;
  cout<<"Enter fesc:"<<endl;
  cin>>fesc;
  cout<<"Enter nion:"<<endl;
  cin>>nion;
  cout<<"Enter fx:"<<endl;
  cin>>fx;
  cout<<"Enter nlya:"<<endl;
  cin>>nlya;
  cout<<"Enter popflag:"<<endl;
  cin>>popflag;
  cout<<"Enter xrayflag:"<<endl;
  cin>>xrayflag;
  cout<<"Enter lyaxrayflag:"<<endl;
  cin>>lyaxrayflag;
  cout<<"Enter sourceflag:"<<endl;
  cin>>sourceflag;

  a.initAstrophysics(fstar,fesc,nion,fx,nlya,popflag,xrayflag,lyaxrayflag,sourceflag,-1.0);

  lyf.initLyman(&c,&a);



////////////////////////////////////////////////////////////////////////
// Calculate history for single model
///////////////////////////////////////////////////////////////////////
     string data_dir;
     string file_like;
     string files;
     double zeta0, zeta1,z0,deltaz;
     double zeta2, zeta3;
     double like, z;
     vector<double> param;
     Ionization Ion(&c,&a);
     double zstep;
     double Nion0(1.0e51);
     double p1, p2, p3, p4;
     double clump, Nion, xi, zeta, mfp, gamma;
     stringstream ss;

     data_dir="./";
     /////////////////////////////////////////////////
     //set up range of parameters
       param_flag=0;
       p1=45.0;
       p2=0.0;
       p3=0.0;
       p4=0.0;

       Ion.setFlags(ndot_flag,cmb_flag,lls_flag,nthread);
       Ion.setGammaFlag(gamma_flag);
       lyf.setLLSFlag(lls_flag);

 

     //core code to get likelihood and output value
     param.clear();
     param.push_back(param_flag);
     param.push_back(p1);
     param.push_back(p2);
     param.push_back(p3);
     param.push_back(p4);
     Ion.setZetaParam(0.0,param,1);

       files=data_dir+"zeta_history";
       ss<<files<<"_"<<p1<<".dat";
       ss>>files;
       fout.open(files.c_str());

       fout<<"#"<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<endl;
       fout<<"# tau="<<Ion.getTauCMB()<<endl;
       fout<<"# z Nion  xi  gamma  zeta  mfp  clump"<<endl;

       z=30.0;
       zstep=0.25;
       while(z>2.5){
	 Nion=Ion.getNion(z)*MPC*MPC*MPC;
	 xi=Ion.getXI(z);
	 gamma=lyf.getGammaFromNdot(z,Nion);
	 zeta=Ion.getZeta(z);
	 mfp=lyf.mfpFromGamma(z,gamma);
	 clump=Ion.getClumping(z,xi);

	 cout<<z<<"\t"<<Nion/Nion0<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<mfp<<"\t"<<clump<<endl;
	 fout<<z<<"\t"<<Nion/Nion0<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<mfp<<"\t"<<clump<<endl;	
	 z-=zstep;
       }

   fout.close();
 



///////////////////////////////////////////////////////////////////////
  return 0;
}




