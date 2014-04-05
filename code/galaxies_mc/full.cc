
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
#include "reionization.h"
#include "haloDensity.h"
#include "spline.h"

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

  // Handle arguments 
  while ((argc>1) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'z':
      zin_arg=atof(&argv[1][2]);
      break;
    case 't':
      tin=atoi(&argv[1][2]);
      break;
    case 'x':
      xin=atoi(&argv[1][2]);
      break;
    case 'l':
      lyaxray_in=atoi(&argv[1][2]);
      break;
    case 'p':
      popflag_in=atoi(&argv[1][2]);
      break;
    default:
      cerr << "Bad option" <<argv[1] <<"\n";
    }
    --argc;
    ++argv;
  }

  
  s8=0.77;
  omb=0.044;
  om0=0.26;
  lam0=0.74;
  h=0.72;
  n=0.95;

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
  
  //////////////////////////////////////////////////////////////
  // Test photo ionisation cross section
  /////////////////////////////////////////////////////////////
  /*
  double nu;
  double sigma;
  double x,z;
  double sigmaHI, sigmaHeI, sigmaHeII;

  file="pi_xsection.dat";
  fout.open(file);
  fout<<"# nu/nuLL, sigma, g1F, HI, HeI, HeII"<<endl;
  nu=1.00001*NULYMANLIMIT/4.0;
  while(nu<=2000.0*NULYMANLIMIT){
    sigma=a.xsectionPhotoIonise(nu,1.0);
    sigmaHI=a.xsectionPhotoIoniseGen(nu,1,1);
    sigmaHeI=a.xsectionPhotoIoniseGen(nu,2,1);
    sigmaHeII=a.xsectionPhotoIoniseGen(nu,2,2);
    x=NULYMANLIMIT/nu;
    z=sqrt(NULYMANLIMIT/(nu-NULYMANLIMIT));
    fout<<nu/NULYMANLIMIT<<"\t"<<sigma<<"\t"<<a.gaunt1F(x,z)<<"\t"<<sigmaHI<<"\t"<<sigmaHeI<<"\t"<<sigmaHeII<<endl;
    nu*=sqrt(sqrt(2.0));
  }
  fout.close();
  */
  /////////////////////////////////////////////////////////////

  /*
  //////////////////////////////////////////////////////////////
  // Test photoionisation mean free path
  /////////////////////////////////////////////////////////////
  double nu,energy;
  double mfp;
  double z,xfrac;

  file="pi_mfp.dat";
  fout.open(file);
  fout<<"# nu/nuLL, sigma, g1F"<<endl;
  nu=1.00001*NULYMANLIMIT;
  z=10.0;
  xfrac=0.1;
  while(nu<=2000.0*NULYMANLIMIT){
    energy=nu/NULYMANLIMIT*RYDBERG;
    mfp=a.mfpPhotoIonise(energy,z,xfrac);
    fout<<nu/NULYMANLIMIT<<"\t"<<mfp<<"\t"<<LIGHTSPEED/c.hubbleZ(z)/UNH/KPC/1.0e3<<endl;
    nu*=sqrt(sqrt(2.0));
  }
  fout.close();
  */
  /////////////////////////////////////////////////////////////
 
  //////////////////////////////////////////////////////////////
  // Test Starburst routines - luminosity and flux
  /////////////////////////////////////////////////////////////
  /*
  double nu,energy;
  double mfp;
  double z,xfrac;
  double r;

  file="starburst_lnu.dat";
  fout.open(file);
  fout<<"# nu/nuLL, sigma, g1F"<<endl;
  nu=1.00001*NULYMANLIMIT;
  z=10.0;
  r=0.1;
  xfrac=0.1;
  while(nu<=2000.0*NULYMANLIMIT){
    energy=nu/NULYMANLIMIT*RYDBERG;
    fout<<nu/NULYMANLIMIT<<"\t"<<a.lumNuStarBurst(nu)<<"\t"<<a.jNuStarBurst(nu,0.1,z)<<"\t"<<a.jNuStarBurst(nu,1.0,z)<<endl;
    nu*=sqrt(sqrt(2.0));
  }
  fout.close();
  */
  ///////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Test Starburst routines - heating
  /////////////////////////////////////////////////////////////
  /*
  double nu,energy;
  double mfp;
  double z,xfrac;
  double r;

  file="starburst_heat.dat";
  fout.open(file);
   nu=1.00001*NULYMANLIMIT;
  z=10.0;
  r=1.0e-3;
  xfrac=0.1;
  fout<<"# SFR="<<"\t"<<a.globalSFR(z)<<endl;
  fout<<"# r, heating rate"<<endl;
  while(r<=1000.0){
    energy=nu/NULYMANLIMIT*RYDBERG;
    fout<<r<<"\t"<<a.heatingStarBurst(r,z)<<endl;
    r*=sqrt(sqrt(2.0));
  }
  fout.close();


  file="starburst_heat_nu.dat";
  fout.open(file);
  nu=1.00001*NULYMANLIMIT;
  z=10.0;
  r=1.0e-3;
  xfrac=0.1;
  fout<<"# SFR="<<"\t"<<a.globalSFR(z)<<endl;
  fout<<"# nu/nuLL, r=0.1Mpc, r=1.0Mpc"<<endl;
  while(nu<=2000.0*NULYMANLIMIT){
    energy=nu/NULYMANLIMIT*RYDBERG;
    fout<<nu/NULYMANLIMIT<<"\t"<<a.heatingNuStarBurst(nu,0.1,z)<<"\t"<<a.heatingNuStarBurst(nu,1.0,z)<<endl;
    nu*=sqrt(sqrt(2.0));
  }
  fout.close();
  */ 
  //////////////////////////////////////////////////////////////
  // Test Shull and van Steenberg fractions
  /////////////////////////////////////////////////////////////
  /*
  double x,f1,f2,f3,f4,f5;
  double xstep;

  file="fractionSVS.dat";
  fout.open(file);
  fout<<"# r, heating rate"<<endl;
  x=1.0e-4;

  xstep=exp(log(1.0/1.0e-1)/10.0);

  while(x<=1.1){
    f1=a.fracPrimaryElectron(x,1);
    f2=a.fracPrimaryElectron(x,2);
    f3=a.fracPrimaryElectron(x,3);
    f4=a.fracPrimaryElectron(x,4);
    f5=a.fracPrimaryElectron(x,5);
    fout<<x<<"\t"<<f1<<"\t"<<f2<<"\t"<<f3<<"\t"<<f4<<"\t"<<f5<<endl;
    x*=xstep;
  }
  fout.close();
*/

  //////////////////////////////////////////////////////////////
  // Test star formation 
  /////////////////////////////////////////////////////////////
  /*
  double z,sfr;

  file="starformationrate.dat";
  fout.open(file);
  fout<<"# z, sfr"<<endl;
  z=0.01;
  while(z<40.0){
    sfr=a.globalSFR(z,1);
    fout<<z<<"\t"<<sfr<<"\t"<<a.globalSFR(z,2)<<endl;
    cout<<z<<"\t"<<sfr<<endl;
    z+=1.0;
  }
  fout.close();
 
  */

  //////////////////////////////////////////////////////////////
  // Test Starburst routines - heating - in Kelvin
  /////////////////////////////////////////////////////////////
  /*
  double nu,energy;
  double mfp;
  double z,xfrac;
  double r;
  double CT,H;

  file="starburst_heat_K.dat";
  fout.open(file);
  nu=1.00001*NULYMANLIMIT;
  z=10.0;
  r=1.0e-3;
  xfrac=0.1;
  fout<<"# SFR="<<"\t"<<a.globalSFR(z)<<endl;
  fout<<"# r, heating rate, DT"<<endl;
  CT=1.5*BOLTZK*c.nh(z);   //heat capacity
  H=c.hubbleZ(z)*UNH*c.getH();
  cout<<BOLTZK<<"\t"<<CT<<"\t"<<H<<endl;

  while(r<=100.0){
    energy=nu/NULYMANLIMIT*RYDBERG;
    fout<<r<<"\t"<<a.heatingStarBurst(r,z)<<"\t"<<a.heatingStarBurst(r,z)/CT/H<<endl;
    r*=sqrt(sqrt(2.0));
  }
  fout.close();

  file="starburst_heat_nu_K.dat";
  fout.open(file);
  fout<<"# SFR="<<"\t"<<a.globalSFR(z)<<endl;
  fout<<"# nu/nuLL, r=0.1Mpc, r=1.0Mpc"<<endl;
  nu=1.00001*NULYMANLIMIT;
  z=10.0;
  r=1.0e-3;
  xfrac=0.1;
  while(nu<=2000.0*NULYMANLIMIT){
    energy=nu/NULYMANLIMIT*RYDBERG;
    fout<<nu/NULYMANLIMIT<<"\t"<<a.heatingNuStarBurst(nu,0.1,z)/CT/H<<"\t"<<a.heatingNuStarBurst(nu,1.0,z)/CT/H<<endl;
    nu*=sqrt(sqrt(2.0));
  }
  fout.close();
 
  */
  //////////////////////////////////////////////////////////////
  // Calculate the global temperature evolution
  /////////////////////////////////////////////////////////////
    /*
    double Tgas,z;

    z=6.0;
    file="Thistory.dat";
    fout.open(file);
    while(z<=36.0){
      Tgas=a.globalTGas(z);
      fout<<z<<"\t"<<Tgas<<endl;
      cout<<z<<"\t"<<Tgas<<endl;
      z+=1.0;
    }
    fout.close();
    */

  //////////////////////////////////////////////////////////////
  // Calculate the global evolution of T & x_e
  /////////////////////////////////////////////////////////////
  /*
    double Tgas,z,xe;
    double *y;

    y=dvector(1,2);

    
    z=6.0;
    file="Thistory2.dat";
    fout.open(file);
    while(z<=36.0){
      a.globalHistory(z,y);
      Tgas=y[1];
      xe=y[2];
      cout<<z<<"\t"<<Tgas<<"\t"<<xe<<endl;
      fout<<z<<"\t"<<Tgas<<"\t"<<xe<<endl;
      z+=1.0;
    }
    fout.close();

    free_dvector(y,1,2);

    double *zs,*xes,trash;
    int zvar(30);

    //get tau value
    zs=dvector(1,zvar);
    xes=dvector(1,zvar);

    file="Thistory2.dat";
    fin.open(file);
    for(i=1;i<=zvar;i++){
      fin>>zs[i]>>trash>>xes[i];
      cout<<zs[i]<<"\t"<<xes[i]<<endl;
    }
    fin.close();

    a.splineXFree(z,zs,xes,zvar,0);

    file="Thistory2spline.dat";
    fout.open(file);
    while(z<=36.0){
      xe=a.getXFree(z);
      cout<<z<<"\t"<<xe<<endl;
      fout<<z<<"\t"<<xe<<endl;
      z+=0.1;
    }
    fout.close();
      
  */
  //////////////////////////////////////////////////////////////
  // Calculate the Optical depth
  /////////////////////////////////////////////////////////////
  /*
  int zvar(0);

  double Tgas,z,xe;
  double *y;
  
  y=dvector(1,2);
  
  z=6.0;
  a.globalHistory(z,y);
  Tgas=y[1];
  xe=y[2];
  cout<<z<<"\t"<<Tgas<<"\t"<<xe<<endl;
  z+=1.0;

  free_dvector(y,1,2);
 
  file="Thistory_save.dat";
  cout<<a.getTauCMB(file)<<endl;
  */
  //////////////////////////////////////////////////////////////
  // Calculate the Clumping Factor, MEHR (2000) - P_V
  /////////////////////////////////////////////////////////////
  /*
  double z;
  double clump;
  double delta;

  file="PVplot_z2.dat";
  fout.open(file);
  z=2.0;  
  a.initPVParamMEHR(z);
  delta=1.0e-2;
  while(delta<100.0){
    fout<<delta<<"\t"<<getPVMEHR(delta)<<"\t"<<getDeltaPVMEHR(delta)<<"\t"<<getDelta2PVMEHR(delta)<<endl;
    delta*=1.2;
  }
  fout.close();

  file="PVplot_z3.dat";
  fout.open(file);
  z=3.0;  
  a.initPVParamMEHR(z);
  delta=1.0e-2;
  while(delta<100.0){
    fout<<delta<<"\t"<<getPVMEHR(delta)<<"\t"<<getDeltaPVMEHR(delta)<<"\t"<<getDelta2PVMEHR(delta)<<endl;
    delta*=1.2;
  }
  fout.close();

  file="PVplot_z4.dat";
  fout.open(file);
  z=4.0;  
  a.initPVParamMEHR(z);
  delta=1.0e-2;
  while(delta<100.0){
    fout<<delta<<"\t"<<getPVMEHR(delta)<<"\t"<<getDeltaPVMEHR(delta)<<"\t"<<getDelta2PVMEHR(delta)<<endl;
    delta*=1.2;
  }
  fout.close();

  file="PVplot_z6.dat";
  fout.open(file);
  z=6.0;  
  a.initPVParamMEHR(z);
  delta=1.0e-2;
  while(delta<100.0){
    fout<<delta<<"\t"<<getPVMEHR(delta)<<"\t"<<getDeltaPVMEHR(delta)<<"\t"<<getDelta2PVMEHR(delta)<<endl;
    delta*=1.2;
  }
  fout.close();
  */
  //////////////////////////////////////////////////////////////
  // Calculate the Clumping Factor, MEHR (2000) - R/R_u
  /////////////////////////////////////////////////////////////
  /*
  double Deltai,z;
  double DeltaMin(1.0e-4);
  double tol(1.0e-4);
  double Fv, Fm, clump;
  file="clumpMEHR_z2.dat";
  z=2.0;
  fout.open(file);
  Deltai=1.0e-2;
  while(Deltai<=1.0e3){
    a.initPVParamMEHR(z);
    Fv=qromb(getPVMEHR,DeltaMin,Deltai,tol);
    Fm=qromb(getDeltaPVMEHR,DeltaMin,Deltai,tol);
    clump=qromb(getDelta2PVMEHR,DeltaMin,Deltai,tol);
    fout<<Deltai<<"\t"<<Fv<<"\t"<<Fm<<"\t"<<clump<<endl;
    Deltai*=1.2;
  }
  fout.close();

  file="clumpMEHR_z3.dat";
  z=3.0;
  fout.open(file);
  Deltai=1.0e-2;
  while(Deltai<=1.0e3){
    a.initPVParamMEHR(z);
    Fv=qromb(getPVMEHR,DeltaMin,Deltai,tol);
    Fm=qromb(getDeltaPVMEHR,DeltaMin,Deltai,tol);
    clump=qromb(getDelta2PVMEHR,DeltaMin,Deltai,tol);
    fout<<Deltai<<"\t"<<Fv<<"\t"<<Fm<<"\t"<<clump<<endl;
    Deltai*=1.2;
  }
  fout.close();

  file="clumpMEHR_z4.dat";
  z=4.0;
  fout.open(file);
  Deltai=1.0e-2;
  while(Deltai<=1.0e3){
    a.initPVParamMEHR(z);
    Fv=qromb(getPVMEHR,DeltaMin,Deltai,tol);
    Fm=qromb(getDeltaPVMEHR,DeltaMin,Deltai,tol);
    clump=qromb(getDelta2PVMEHR,DeltaMin,Deltai,tol);
    fout<<Deltai<<"\t"<<Fv<<"\t"<<Fm<<"\t"<<clump<<endl;
    Deltai*=1.2;
  }
  fout.close();

  file="clumpMEHR_z6.dat";
  z=6.0;
  fout.open(file);
  Deltai=1.0e-2;
  while(Deltai<=1.0e3){
    a.initPVParamMEHR(z);
    Fv=qromb(getPVMEHR,DeltaMin,Deltai,tol);
    Fm=qromb(getDeltaPVMEHR,DeltaMin,Deltai,tol);
    clump=qromb(getDelta2PVMEHR,DeltaMin,Deltai,tol);
    fout<<Deltai<<"\t"<<Fv<<"\t"<<Fm<<"\t"<<clump<<endl;
    Deltai*=1.2;
  }
  fout.close();

  */

  //////////////////////////////////////////////////////////////
  // Calculate 21cm collisional rates
  /////////////////////////////////////////////////////////////
  /*
  double Tk;
  double  kappaEH,kappaHH;

  file="kappa_coll.dat";
  fout.open(file);

  Tk=1.0e0;
  while(Tk<1.0e4){
    kappaEH=tocm.kappaEH(Tk);
    kappaHH=tocm.kappaHH(Tk);
    fout<<Tk<<"\t"<<kappaEH<<"\t"<<kappaHH<<endl;
    Tk*=1.5;
  }
  fout.close();
 
  */
  //////////////////////////////////////////////////////////////
  // Calculate 21cm collisional coupling: Nasser (2005) Fig 2.
  /////////////////////////////////////////////////////////////
  /*
  double Tk, xi,z;
  double xc,tb,ts;

  file="xcoll_z15_T4.dat";
  fout.open(file);

  z=15.0;
  Tk=1.0e4;
  //  while(Tk>200.0){
    xi=1.0;   
    while(xi>-0.01){
      xc=tocm.getXColl(z,Tk,xi);
      tb=tocm.tBrightGen(z,Tk,xi,0.0);
      ts=tocm.tSpinGen(z,Tk,xi,0.0);
      fout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      cout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      xi-=0.025;
    }
    //  Tk-=250.0;
    //}
  fout.close();

  file="xcoll_z20_T4.dat";
  fout.open(file);

  z=20.0;
  Tk=1.0e4;
  //  while(Tk>200.0){
    xi=1.0;   
    while(xi>-0.01){
      xc=tocm.getXColl(z,Tk,xi);
      tb=tocm.tBrightGen(z,Tk,xi,0.0);
      ts=tocm.tSpinGen(z,Tk,xi,0.0);
      fout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      cout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      xi-=0.025;
    }
    //  Tk-=250.0;
    //}
  fout.close();

  file="xcoll_z15_T3.dat";
  fout.open(file);

  z=15.0;
  Tk=1.0e3;
  //  while(Tk>200.0){
    xi=1.0;   
    while(xi>-0.01){
      xc=tocm.getXColl(z,Tk,xi);
      tb=tocm.tBrightGen(z,Tk,xi,0.0);
      ts=tocm.tSpinGen(z,Tk,xi,0.0);
      fout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      cout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      xi-=0.025;
    }
    //  Tk-=250.0;
    //}
  fout.close();

  file="xcoll_z20_T3.dat";
  fout.open(file);

  z=20.0;
  Tk=1.0e3;
  //  while(Tk>200.0){
    xi=1.0;   
    while(xi>-0.01){
      if(xi<0.0) xi=0.0;
      xc=tocm.getXColl(z,Tk,xi);
      tb=tocm.tBrightGen(z,Tk,xi,0.0);
      ts=tocm.tSpinGen(z,Tk,xi,0.0);
      fout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      cout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      xi-=0.025;
    }
    //  Tk-=250.0;
    //}
  fout.close();
  */

  //////////////////////////////////////////////////////////////
  // Calculate 21cm collisional coupling: Nasser (2005) Fig 1.
  /////////////////////////////////////////////////////////////
  /*
  double Tk, xi,z;
  double xc,tb,ts;

  file="xcollmap_z15_T4.dat";
  fout.open(file);

  z=15.0;
  Tk=1.0e4;
  while(Tk>200.0){
    xi=1.0;   
    while(xi>-0.01){
      xc=tocm.getXColl(z,Tk,xi);
      tb=tocm.tBrightGen(z,Tk,xi,0.0);
      ts=tocm.tSpinGen(z,Tk,xi,0.0);
      fout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      cout<<Tk<<"\t"<<xi<<"\t"<<xc<<"\t"<<tb<<"\t"<<ts<<endl;
      xi-=0.025;
    }
    Tk-=250.0;
  }
  fout.close();

  */
  //////////////////////////////////////////////////////////////
  // Test Emissivity code
  /////////////////////////////////////////////////////////////
  /*
  double nu,z;

  file="emission.dat";
  fout.open(file);

  z=10.0;
  nu=3.0/4.0*NULYMANLIMIT;
  while(nu<=NULYMANLIMIT){
    fout<<nu<<"\t"<<tocm.sourceEmission(nu,z)<<"\t"<<a.globalSFR(z)<<endl;
    nu*=1.01;
  }
  fout.close();
  */
  //////////////////////////////////////////////////////////////
  // Test Lya Flux calculation
  /////////////////////////////////////////////////////////////
  /*
  double z, lyaflux;
  double tb,ts,xa,xc;
  double tk,xi,xe,xitot;
  double *result;

  result=dvector(1,3);

  file="lyaflux.dat";
  fout.open(file);
  z=35.0;
  while(z>5.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];  
    xe=result[3];  
    lyaflux=tocm.lyaFlux(z);
    xc=tocm.getXColl(z,tk,xi);
    xa=tocm.getXAlpha(z,tk,lyaflux);
    ts=tocm.tSpinGen(z,tk,xi,lyaflux);
    tb=(1.0-xi)*tocm.tBrightGen(z,tk,xe,lyaflux);
    xitot=xi*xi+(1.0-xi)*xe;
    
    cout<<z<<"\t"<<tk<<"\t"<<xi<<"\t"<<lyaflux<<"\t"<<xc<<"\t"<<xa<<"\t"<<ts<<"\t"<<tb*1.0e3<<"\t"<<xe<<"\t"<<xitot<<endl;
    fout<<z<<"\t"<<tk<<"\t"<<xi<<"\t"<<lyaflux<<"\t"<<xc<<"\t"<<xa<<"\t"<<ts<<"\t"<<tb*1.0e3<<"\t"<<xe<<"\t"<<xitot<<endl;
    z-=0.5;
  }
  fout.close();

  free_dvector(result,1,3);
  
  */
  //////////////////////////////////////////////////////////////
  // Calculate the global temperature evolution - three variable
  /////////////////////////////////////////////////////////////
  /*
    double Tgas,z;
    double *result;
    result=dvector(1,2);

    z=6.0;
    //file="Thistory.dat";
    //fout.open(file);
    //while(z>=5){
      
    a.globalHistory(z,result);
    //   fout<<z<<"\t"<<Tgas<<endl;
    //cout<<z<<"\t"<<Tgas<<endl;
      //     z-=1.0;
      //}
      //fout.close();
  
    free_dvector(result,1,2);
  */
  //////////////////////////////////////////////////////////////
  // Test RECFAST calling code
  /////////////////////////////////////////////////////////////
  /*
  double z, *result;

  result=dvector(1,3);

  z=40.0;
  // a.callRECFAST();
  //a.useRECFAST(z,result);
  a.globalHistory(20.0,result);
  a.getTIGM(20.0,result);

  free_dvector(result,1,3);
  */
  //////////////////////////////////////////////////////////////
  // Test Lya fluctuation code - window function
  /////////////////////////////////////////////////////////////
  /*
  double z,k,wk;

  double *result;
  double m,dsdm;
  double zuse,dJdz,tk,lyaflux,xa;
  double xanorm;
  
  z=20.0;
  result=dvector(1,3);
  a.getTIGM(z,result);
  tk=result[1];
  lyaflux=tocm.lyaFlux(z);

  xa = tocm.getXAlpha(z,tk,lyaflux);
  xanorm=1.0/xa;
  cout<<tk<<"\t"<<xa<<endl;
 
  file="lyawindow.dat";
  fout.open(file);

  zuse=z;
  k=0.001;
  setWindowK(z,zuse,k,n,tk,&c,&a,&tocm,1);
  while(zuse<tocm.lymanZMax(z,2.0)+0.12){
    dJdz=setWindowLyaK(z,zuse,k,n,tk,&c,&a,&tocm,0)*xanorm;
    cout<<zuse<<"\t"<<dJdz<<endl;
    fout<<zuse<<"\t"<<dJdz<<endl;
    zuse+=0.001;
  }
  fout.close();
  
  file="lyaprofile.dat";
  fout.open(file);

  zuse=z;
  while(zuse<tocm.lymanZMax(z,2.0)+0.12){
    n=2.0;
    dJdz=0.0;
    while((zuse<=tocm.lymanZMax(z,n))&&(int(n)<=lynmax)){
      dJdz+=setDJalphaDz(n,z,zuse,&c,&a,&tocm,1)*tocm.lyaRecycle(n)*xanorm;
      n+=1.0;
    }
    cout<<zuse<<"\t"<<dJdz<<endl;
    fout<<zuse<<"\t"<<dJdz<<endl;
    zuse+=0.001;
  }
  fout.close();

 
  file="lyafluctuate.dat";
  fout.open(file);
  sigm(m,dsdm,&c,2);
  z=20.0;
  k=1.0e-3;
  while(k<2.0e3){
    wk=tocm.getWindowLya(z,k);
    //    cout<<k<<"\t"<<wk<<endl;
    fout<<k<<"\t"<<wk<<endl;
    k*=2.0;
  }
  fout.close();
  
  */
  //////////////////////////////////////////////////////////////
  // Test Lya fluctuation coefficients
  /////////////////////////////////////////////////////////////
  /*
  double z;
  double *result,*beta;
  double tk,lyaflux,xi,xe;
  double xc,xa,ts,tb,xitot;

  result=dvector(1,3);
  beta=dvector(1,5);

  file="beta_history.dat";
  fout.open(file);

  z=40.0;
  while(z>6.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    lyaflux=tocm.lyaFlux(z);
    xc=tocm.getXColl(z,tk,xi);
    xa=tocm.getXAlpha(z,tk,lyaflux);
    ts=tocm.tSpinGen(z,tk,xi,lyaflux);
    tb=(1.0-xi)*tocm.tBrightGen(z,tk,xe,lyaflux);
    xitot=xi*xi+(1.0-xi)*xe;
    tocm.getBeta(z,tk,xe,lyaflux,beta);
    z-=0.5;
    cout<<z<<"\t"<<tk<<"\t"<<xe<<"\t"<<lyaflux<<"\t"<<beta[1]<<"\t"<<beta[2]<<"\t"<<beta[3]<<"\t"<<beta[4]<<"\t"<<tb<<endl;
    fout<<z<<"\t"<<tk<<"\t"<<xe<<"\t"<<lyaflux<<"\t"<<beta[1]<<"\t"<<beta[2]<<"\t"<<beta[3]<<"\t"<<beta[4]<<"\t"<<tb<<endl;
  }

  fout.close();

  free_dvector(result,1,3);
  free_dvector(beta,1,5);
  */
  //////////////////////////////////////////////////////////////
  // Test 21cm coupling coefficients
  /////////////////////////////////////////////////////////////
  /*
  double z;
  double *result,*beta;
  double tk,lyaflux,xi,xe;
  double xc,xa,ts,tb,xitot;
  double xcEH,xcHH;

  result=dvector(1,3);
  beta=dvector(1,5);

  file="xc_history.dat";
  fout.open(file);

  z=40.0;
  while(z>6.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    lyaflux=tocm.lyaFlux(z);
    xc=tocm.getXColl(z,tk,xi);
    xcEH=tocm.getXCollEH(z,tk,xi);
    xcHH=tocm.getXCollHH(z,tk,xi);
    xa=tocm.getXAlpha(z,tk,xe,lyaflux);
    ts=tocm.tSpinGen(z,tk,xi,lyaflux);
    tb=(1.0-xi)*tocm.tBrightGen(z,tk,xe,lyaflux);
    xitot=xi*xi+(1.0-xi)*xe;
    tocm.getBeta(z,tk,xe,lyaflux,beta);
    z-=0.5;
    cout<<z<<"\t"<<tk<<"\t"<<xe<<"\t"<<lyaflux<<"\t"<<xa<<"\t"<<xc<<"\t"<<xcEH<<"\t"<<xcHH<<"\t"<<tb<<endl;
    fout<<z<<"\t"<<tk<<"\t"<<xe<<"\t"<<lyaflux<<"\t"<<xa<<"\t"<<xc<<"\t"<<xcEH<<"\t"<<xcHH<<"\t"<<tb<<endl;
  }

  fout.close();

  free_dvector(result,1,3);
  free_dvector(beta,1,5);
 
  */
 //////////////////////////////////////////////////////////////
  // Test Lya fluctuation code - power spectrum
  /////////////////////////////////////////////////////////////
  /*
  double z,k,pkfull;

  double *result,*store,*beta;
  double m,dsdm;
  double zuse,dJdz,tk,lyaflux,xa;
  double rT,rF,xi,xe;

  extern double xanorm;
  
  z=14.0;
  result=dvector(1,3);
  store=dvector(1,5);
  beta=dvector(1,5);

  a.getTIGM(z,result);
  tk=result[1];
  xi=result[2];
  xe=result[3];
  lyaflux=tocm.lyaFlux(z);

  rF=tocm.getCutoffFilter(z);
  rT=tocm.getCutoffThermal(z,tk);

  xa = tocm.getXAlpha(z,tk,lyaflux);

  tocm.getBeta(z,tk,xi,lyaflux,beta);  //xi or xe?
  cout<<beta[1]<<"\t"<<beta[2]<<"\t"<<beta[3]<<"\t"<<beta[4]<<endl;


  // xanorm=xa;
  //lyaflux=tocm.lyaFlux(z);
  //xa = tocm.getXAlpha(z,tk,lyaflux);
  cout<<tk<<"\t"<<xa<<"\t"<<rT<<"\t"<<rF<<endl;
  
  file="lyapower.dat";
  fout.open(file);
  sigm(m,dsdm,&c,2);

  k=1.0e-3;
  while(k<2.0e3){
    pkfull=tocm.powerSpectrumLya(z,k,store);
    //    cout<<k<<"\t"<<wk<<endl;
    fout<<k<<"\t"<<store[1]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    k*=2.0;
  }
  fout.close();
  
  free_dvector(beta,1,5);
  free_dvector(result,1,3);
  free_dvector(store,1,5);
*/
 //////////////////////////////////////////////////////////////
  // Characterise Lya excess power as function of z
  /////////////////////////////////////////////////////////////
  /*
  
 
  double z,k,pkfull;

  double *result,*store,*beta;
  double m,dsdm;
  double zuse,dJdz,tk,lyaflux,xa;
  double rT,rF,xi,xe;
  double gammaA,betabar;

  extern double xanorm;
  
  result=dvector(1,3);
  store=dvector(1,5);
  beta=dvector(1,5);


  
  file="lyapowerZ.dat";
  fout.open(file);

  z=25.0;
  k=1.0e-0;

  fout<<"#"<<k<<endl;
  while(z>6.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    lyaflux=tocm.lyaFlux(z);
    
    rF=tocm.getCutoffFilter(z);
    rT=tocm.getCutoffThermal(z,tk);
    
    xa = tocm.getXAlpha(z,tk,lyaflux);
    
    tocm.getBetaTb(z,tk,xi,lyaflux,beta);  //xi or xe?
    cout<<beta[1]<<"\t"<<beta[2]<<"\t"<<beta[3]<<"\t"<<beta[4]<<"\t"<<beta[5]<<endl;
    gammaA=tocm.getGammaA(z,tk);
    betabar=beta[1]+beta[4]*gammaA;

    cout<<tk<<"\t"<<xa<<"\t"<<rT<<"\t"<<rF<<endl;
    
    pkfull=tocm.powerSpectrumLya(z,k,store);
    cout<<z<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<"\t"<<betabar<<"\t"<<betabar+store[1]*beta[3]<<endl;
    fout<<z<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<"\t"<<betabar<<"\t"<<betabar+store[1]*beta[3]<<endl;
    z-=0.5;
  }
  fout.close();
  
  free_dvector(beta,1,5);
  free_dvector(result,1,3);
  free_dvector(store,1,5);
  */
 //////////////////////////////////////////////////////////////
  // Investigate heating fluctuations
  /////////////////////////////////////////////////////////////
  /*
 
  double Q,ntot,CT,dzdt,heat;
  double tk,xi,xe,*result;
  double z,zz;
  double Mpc(KPC*1.0e3);
  double QC, QX;

  result=dvector(1,3);

  file="xrayQ.dat";
  fout.open(file);

  z=1000.0;
  while(z>5.0){

    zz=z+1.0;
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    ntot=c.nh(z)*(1.0+xe+FHE);  //n=n_H+n_e+n_He
    CT=1.5*BOLTZK*ntot;  //specific heat capacity in ergs K^{-1} cm^{-3}
    //  CT=2.5*BOLTZK*ntot;  //specific heat capacity in ergs K^{-1} cm^{-3}
    dzdt=zz*c.hubbleZ(z)*UNH*c.getH();
    
    heat=a.lumTotStarBurst(a.globalSFR(z))/Mpc/Mpc/Mpc*a.fracPrimaryElectron(xe,1); //xrays
    heat/=dzdt*CT; 
    
    Q=heat/tk;  //Coefficient on RHS of delta_T evolution equation
    
    QC=tocm.QCompton(z);
    QX=tocm.QXray(z);

    cout<<z<<"\t"<<tk<<"\t"<<heat<<"\t"<<Q<<"\t"<<QC<<"\t"<<QX<<endl;
    fout<<z<<"\t"<<tk<<"\t"<<heat<<"\t"<<Q<<"\t"<<QC<<"\t"<<QX<<endl;    

    z-=0.5;
  }

  fout.close();
  free_dvector(result,1,3);

  */
  //////////////////////////////////////////////////////////////
  // Test mean x-ray code - total heating rate
  /////////////////////////////////////////////////////////////
  // Compare heating rate from flux with heating rate assuming
  // all energy is deposited into IGM 
  //
  /*

  double flux,heat;
  double fluxC,heatC;
  double z,zz;

  file="xrayflux.dat";
  fout.open(file);

  z=40.0;
  while(z>5.0){
    zz=z+1.0;
    fluxC=a.lumTotStarBurst(a.globalSFR(z))/zz/zz/zz;
    fluxC/=1.0e3/RYDBERG*RYDBERG_CGS;
    fluxC/=4.0*PI*MPC*MPC*1.0e2; //flux at distance 10Mpc

    heatC=a.lumTotStarBurst(a.globalSFR(z));
    heatC/=MPC*MPC*MPC;
    //    heatC/=zz*zz*zz*zz*zz;

    flux=tocm.xrayFlux(z);
    heat=tocm.xrayHeating(z);

    cout<<z<<"\t"<<flux<<"\t"<<fluxC<<"\t"<<heat<<"\t"<<heatC<<endl;
    fout<<z<<"\t"<<flux<<"\t"<<fluxC<<"\t"<<heat<<"\t"<<heatC<<endl; 

    //    cout<<RYDBERG/RYDBERG_CGS<<endl;

    z-=0.5;
  }
  
  fout.close();

  //////////////
  double zp,xrayheat;

  file="xrayfluxZ_z10.dat";
  fout.open(file);
  zp=10.0;
  z=zp-1.0e-3;
  while(z>5.0){
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    fout<<z<<"\t"<<zp<<"\t"<<xrayheat<<endl;
    if(fabs(z-zp)<1.0) {
      z-=0.1;
    }else{
      z-=0.5;
    }
  }
  fout.close();

  file="xrayfluxZ_z15.dat";
  fout.open(file);
  zp=15.0;
  z=zp-1.0e-3;
  while(z>5.0){
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    fout<<z<<"\t"<<zp<<"\t"<<xrayheat<<endl;
    if(fabs(z-zp)<1.0) {
      z-=0.1;
    }else{
      z-=0.5;
    }
  }
  fout.close();

  file="xrayfluxZ_z20.dat";
  fout.open(file);
  zp=20.0;
  z=zp-1.0e-3;
  while(z>5.0){
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    fout<<z<<"\t"<<zp<<"\t"<<xrayheat<<endl;
    if(fabs(z-zp)<1.0) {
      z-=0.1;
    }else{
      z-=0.5;
    }
  }
  fout.close();

  file="xrayfluxZ_z25.dat";
  fout.open(file);
  zp=25.0;
  z=zp-1.0e-3;
  while(z>5.0){
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    fout<<z<<"\t"<<zp<<"\t"<<xrayheat<<endl;
    if(fabs(z-zp)<1.0) {
      z-=0.1;
    }else{
      z-=0.5;
    }
  }
  fout.close();

  file="xrayfluxZ_z30.dat";
  fout.open(file);
  zp=30.0;
  z=zp-1.0e-3;
  while(z>5.0){
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    fout<<z<<"\t"<<zp<<"\t"<<xrayheat<<endl;
    if(fabs(z-zp)<1.0) {
      z-=0.1;
    }else{
      z-=0.5;
    }
  }
  fout.close();

  file="xrayfluxZ_z35.dat";
  fout.open(file);
  zp=35.0;
  z=zp-1.0e-3;
  while(z>5.0){
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    fout<<z<<"\t"<<zp<<"\t"<<xrayheat<<endl;
    if(fabs(z-zp)<1.0) {
      z-=0.1;
    }else{
      z-=0.5;
    }
  }
  fout.close();

  file="xrayfluxZ_z40.dat";
  fout.open(file);
  zp=40.0;
  z=zp-1.0e-3;
  while(z>5.0){
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    fout<<z<<"\t"<<zp<<"\t"<<xrayheat<<endl;
    if(fabs(z-zp)<1.0) {
      z-=0.1;
    }else{
      z-=0.5;
    }
  }
  fout.close();

  */
  //////////////////////////////////////////////////////////////
  // Test mean x-ray code - energy deposition with z and r
  /////////////////////////////////////////////////////////////
  /*
  // Compare heating rate from flux with heating rate assuming
  // all energy is deposited into IGM 
  //

  double xrayflux;
  double z,zz;
  double flux1;
  double heat,r;
  double zp,xrayheat;
  double lyaflux,xa;
  double tk,xi,xe,*result;
  result=dvector(1,3);
  file="xrayheatR.dat";
  fout.open(file);

  zp=15.0;
  z=zp-1.0e-6;
  i=0;  //linear close in then log steps later
  while(z>5.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    xrayflux=getDJxrayDz(zp);

    lyaflux=tocm.getDJalphaDz(z,zp);
    xa=tocm.getXAlpha(z,tk,lyaflux);

    r=tocm.relativeConformalR(z,zp);

    cout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<xrayflux<<"\t"<<xrayheat<<"\t"<<lyaflux<<"\t"<<xa<<endl;
    fout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<xrayflux<<"\t"<<xrayheat<<"\t"<<lyaflux<<"\t"<<xa<<endl;

    if(i<30){
      z-=0.0001;
    }else{
      z/=1.001;
    }
    i++;
  }
  
  fout.close();

  free_dvector(result,1,3);
  */
  //////////////////////////////////////////////////////////////
  // Test x-ray optical depth with z and r
  /////////////////////////////////////////////////////////////
  // Compare heating rate from flux with heating rate assuming
  // all energy is deposited into IGM 
  //
  /*
  double z,zz;
  double flux1;
  double heat,r;
  double zp;
  double tk,xi,xe,*result;
  double tau100,tau500,tau1000,tau10000;
  double tau100s,tau500s,tau1000s,tau10000s;
  result=dvector(1,3);
  file="xraytauR.dat";
  fout.open(file);

  zp=25.0;
  z=zp-1.0e-6;
  i=0;  //linear close in then log steps later
  while(z>5.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);

    r=tocm.relativeConformalR(z,zp);
    tau100=tocm.XrayTau(z,zp,100.0);
    tau500=tocm.XrayTau(z,zp,500.0);
    tau1000=tocm.XrayTau(z,zp,1000.0);
    tau10000=tocm.XrayTau(z,zp,10000.0);
    tau100s=tocm.XrayTauSimp(z,zp,100.0);
    tau500s=tocm.XrayTauSimp(z,zp,500.0);
    tau1000s=tocm.XrayTauSimp(z,zp,1000.0);
    tau10000s=tocm.XrayTauSimp(z,zp,10000.0);


    cout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<tau100<<"\t"<<tau500<<"\t"<<tau1000<<"\t"<<tau10000<<"\t"<<tau100s<<"\t"<<tau500s<<"\t"<<tau1000s<<"\t"<<tau10000s<<endl;
    fout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<tau100<<"\t"<<tau500<<"\t"<<tau1000<<"\t"<<tau10000<<"\t"<<tau100s<<"\t"<<tau500s<<"\t"<<tau1000s<<"\t"<<tau10000s<<endl;

    if(i<30){
      z-=0.0001;
    }else{
      z/=1.001;
    }
    i++;
  }
  
  fout.close();

  free_dvector(result,1,3);
  */ 
  //////////////////////////////////////////////////////////////
  // Test x-ray fluctuation code - power spectrum -k  Lambda approx
  /////////////////////////////////////////////////////////////
  /*
  double k,z;
  double *store;
  store=dvector(1,5);

  file="xrayfluctuateL.dat";
  fout.open(file);
  cout<<z<<endl;
  fout<<"#"<<z<<endl;

  z=20.0;
  k=1.0e-3;
  while(k<2.0e3){
    tocm.powerSpectrumXrayL(z,k,store);
    //tocm.powerSpectrumLya(z,k,store);
    fout<<k<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    k*=sqrt(sqrt(2.0));
  }

  fout.close();

  free_dvector(store,1,5);
  */
 
  //////////////////////////////////////////////////////////////
  // Test x-ray fluctuation code - power spectrum - excess Lambda approx
  /////////////////////////////////////////////////////////////
  /*
  double k,z;
  double *store;
  store=dvector(1,5);

  k=10.0;
  file="xrayfluctuateZ.dat";
  fout.open(file);
  cout<<"k="<<k<<endl;
  fout<<"#"<<k<<endl;

  z=25.0;
  while(z>5.0){
    tocm.powerSpectrumXrayL(z,k,store);
    cout<<z<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    fout<<z<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    z-=0.5;
  }

  fout.close();

  free_dvector(store,1,5);
  */

 //////////////////////////////////////////////////////////////
  // Test g_T evolution code
  /////////////////////////////////////////////////////////////
  /*
  double gT;
  double tk,xi,xe,*result;
  double z,zz,k;
  double QC, QX;

  result=dvector(1,3);

  //  file="gT.dat";
  //fout.open(file);

   z=1000.0;
  //while(z>5.0){

  //    zz=z+1.0;
  //  a.getTIGM(z,result);
  //  tk=result[1];
  //  xi=result[2];
  //  xe=result[3];

   z=20.0;
   k=0.001;
   cout<<tocm.getGTHistory(z,k)<<endl;
    // QC=tocm.QCompton(z);
    //QX=tocm.QXray(z);

    //    cout<<z<<"\t"<<tk<<"\t"<<gT<<"\t"<<QC<<"\t"<<QX<<endl;
    //  fout<<z<<"\t"<<tk<<"\t"<<gT<<"\t"<<QC<<"\t"<<QX<<endl;    

    //    z-=1.0;
    //}

    //  fout.close();
  free_dvector(result,1,3);
  */
  //////////////////////////////////////////////////////////////
  // Test x-ray fluctuation code - power spectrum - excess
  /////////////////////////////////////////////////////////////
  /*
  double k,z;
  double *store;
  store=dvector(1,5);

  z=10.0;
  file="xrayfluctuate.dat";
  fout.open(file);
  cout<<"z="<<k<<endl;
  fout<<"#"<<k<<endl;

  k=1.0e-3;
  while(k<2.0e3){
    tocm.powerSpectrumXray(z,k,store);
    cout<<k<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    fout<<k<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    k*=2.0;
  }

  fout.close();

  free_dvector(store,1,5);
  */
  //////////////////////////////////////////////////////////////
  // Test x-ray fluctuation code - power spectrum 
  /////////////////////////////////////////////////////////////
  /*
  double k,z;
  double *store;
  store=dvector(1,5);

  k=10.0;
  file="xrayfluctuateZ.dat";
  fout.open(file);
  cout<<"k="<<k<<endl;
  fout<<"#"<<k<<endl;

  z=25.0;
  while(z>5.0){
    tocm.powerSpectrumXray(z,k,store);
    cout<<z<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    fout<<z<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<endl;
    z-=0.5;
  }

  fout.close();

  free_dvector(store,1,5);
  */

  //////////////////////////////////////////////////////////////
  // Undertstanding W(k) kernel
  /////////////////////////////////////////////////////////////
  /*
  double xrayflux;
  double z,zz;
  double flux1;
  double heat,r;
  double zp,xrayheat;
  double lyaflux,xa;
  double dldz,dxdz,dJdz;
  double rz,zuse,x,oscfun,k;
  double tol(1.0e-4);
  double lEmin(log(XrayEmin)),lEmax(log(XrayEmax));

  extern double zpsave;

  double tk,xi,xe,*result;
  result=dvector(1,3);
  file="wkR.dat";
  fout.open(file);

  z=15.0;
  k=0.1;
  zp=z+1.0e-6;
  i=0;  //linear close in then log steps later
  while(zp<45.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];

    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(10.0,z,&a,1);

    zpsave=zp;
    dldz=qromb(getDLambdaXrayDlE,lEmin,lEmax,tol);

    dJdz=tocm.getDJalphaDz(z,zp);
    dxdz=tocm.getXAlpha(z,tk,dJdz);

    zuse=zp;
    rz=SPEEDOFLIGHT_CGS*(c.confTime(z)-c.confTime(zuse))/H0_CMSMPC/c.getH();
    x=rz*k;
    oscfun = (1.0+tocm.splineBias(zuse))*sin(x)/x;
    oscfun -= 2.0*((3.0/pow(x,3.0)-1.0/x)*sin(x)-3.0*cos(x)/x/x)/3.0;
    oscfun*= c.growthFac(zuse)/c.growthFac(z);

    r=tocm.relativeConformalR(z,zp);

    cout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<dldz<<"\t"<<dxdz<<"\t"<<oscfun<<endl;
    fout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<dldz<<"\t"<<dxdz<<"\t"<<oscfun<<endl;

    if(i<30){
      zp+=0.0001;
    }else{
      zp*=1.001;
    }
    i++;
  }
  
  fout.close();

  free_dvector(result,1,3);
  */

  //////////////////////////////////////////////////////////////
  // Output full power spectrum
  /////////////////////////////////////////////////////////////
  /*
  double k,z;
  double *store;
  int t(2);
  store=dvector(1,7);

  cout<<zin<<"\t"<<tin<<endl;

  z=15.0;
  
  z=zin;
  t=tin;

  file="tbfluctuate.dat";
  fout.open(file);
  cout<<"z="<<z<<endl;
  fout<<"#"<<z<<endl;

  k=1.0e-3;
  while(k<2.0e3){
    tocm.powerSpectrumFull(z,k,store,t);
    cout<<k<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<"\t"<<store[6]<<"\t"<<store[7]<<endl;
    fout<<k<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<"\t"<<store[6]<<"\t"<<store[7]<<endl;
    k*=2.0;
  }

  fout.close();

  free_dvector(store,1,7);
  */
  //////////////////////////////////////////////////////////////
  // Output Kinetic Temperature power spectrum
  /////////////////////////////////////////////////////////////
  /*
  double k,z;
  double *store;
  int t(2);
  store=dvector(1,7);

  cout<<zin<<"\t"<<tin<<endl;

  z=15.0;
  
  z=zin;
  t=tin;

  file="tkfluctuate.dat";
  fout.open(file);
  cout<<"z="<<z<<endl;
  fout<<"#"<<z<<endl;

  k=1.0e-3;
  while(k<2.0e3){
    tocm.powerSpectrumTemp(z,k,store,t);
    cout<<k<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<"\t"<<store[6]<<"\t"<<store[7]<<endl;
    fout<<k<<"\t"<<store[1]<<"\t"<<store[2]<<"\t"<<store[3]<<"\t"<<store[4]<<"\t"<<store[5]<<"\t"<<store[6]<<"\t"<<store[7]<<endl;
    k*=2.0;
  }

  fout.close();

  free_dvector(store,1,7);
 */
  ///////////////////////////////////////////////////////////////////////
  /*
  cout<<"calling heating"<<endl;
  cout<<tocm.xrayHeating(20.0)<<endl;
  cout<<"heating done"<<endl;
  */
  //////////////////////////////////////////////////////////////
  // Check evolution of Naoz and Barkana Imprint
  /////////////////////////////////////////////////////////////
  /*

  double k,z,*zin;
  double **gT;
  int t(2);
  int nz(19);
  zin=dvector(1,nz);
  gT=dmatrix(1,3,1,nz);

  z=10.0;
  for(i=1;i<=nz;i++){
    zin[i]=z;
    cout<<zin[i]<<endl;
    z+=5.0;
  }

  k=0.01;
  tocm.getGTN(zin,k,gT,nz);

  file="scale_gT_ic_k001.dat";
  fout.open(file);
  fout<<"#"<<k<<endl;
  for(i=1;i<=nz;i++){
    cout<<zin[i]<<"\t"<<gT[1][i]<<"\t"<<gT[2][i]<<"\t"<<gT[3][i]<<endl;
    fout<<zin[i]<<"\t"<<gT[1][i]<<"\t"<<gT[2][i]<<"\t"<<gT[3][i]<<endl;
  }
  fout.close();
  */

  //////////////////////////////////////////////////////////////
  // Test code for lya flux from xrays - energy deposition with z and r
  /////////////////////////////////////////////////////////////
  /*
  // Compare heating rate from flux with heating rate assuming
  // all energy is deposited into IGM 
  //

  double xrayflux;
  double z,zz;
  double flux1;
  double heat,r;
  double zp,xrayheat;
  double lyaflux,xa;
  double lyafluxX;
  double tk,xi,xe,*result;
  result=dvector(1,3);
  file="lyafluxR.dat";
  fout.open(file);

  zp=20.0;
  z=zp-1.0e-6;
  i=0;  //linear close in then log steps later
  while(z>5.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];

    //    cout<<z<<"\t"<<tk<<"\t"<<xi<<"\t"<<xe<<endl;
    setDJxrayDzDE(10.0,z,0.0,&c,&a,&tocm,1);
    getDLambdaXrayDzDE(xe,z,&a,1);
    getDJalphaXrayDzDE(xe,z,&a,&c,1);
    setDJalphaDz(z,z,&c,&a,&tocm,1);
    xrayheat=getDLambdaXrayDz(zp);
    xrayheat*=c.nh(z);
    xrayflux=getDJxrayDz(zp);

    tocm.lyaFlux(z);
    lyaflux=getDJalphaDz(zp);
    //lyafluxX=0.0;
    lyafluxX=getDJalphaXrayDz(zp);
    xa=tocm.getXAlpha(z,tk,xe,lyaflux);

    r=tocm.relativeConformalR(z,zp);

    cout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<xrayflux<<"\t"<<lyafluxX<<"\t"<<lyaflux<<"\t"<<xa<<endl;
    fout<<z<<"\t"<<zp<<"\t"<<r<<"\t"<<xrayflux<<"\t"<<lyafluxX<<"\t"<<lyaflux<<"\t"<<xa<<endl;

    if(i<30){
      z-=0.0001;
    }else{
      z/=1.001;
    }
    i++;
  }
  
  fout.close();

  free_dvector(result,1,3);
  */
  ///////////////////////////////////////////////////////////////////////
  // Test GT with x extension
  /////////////////////////////////////////////////////////////////////
  /*
  double k;

  k=0.1;

  tocm.getGTHistory(10.0,k);

  */
  ///////////////////////////////////////////////////////////////////////
  // Calculate z evolution of Temp power spectrum
  /////////////////////////////////////////////////////////////////////
  /*
  double z,k,*zin;
  int nz(21);
  double **store;
  int t;

  t=tin;
  
  zin=dvector(1,nz);
  store=dmatrix(1,8,1,nz);

  z=30.0;
  for(i=nz;i>=1;i--){
    cout<<i<<"\t"<<z<<endl;
    zin[i]=z;
    z-=1.0;
  }

  k=0.1;
  if(t<0){
    tocm.powerSpectrumTempN(zin,k,store,nz);
  }else{
    tocm.powerSpectrumFullN(zin,k,store,nz,t);
  }

  file="tempZ.dat";
  fout.open(file);
  fout<<"#"<<k<<endl;
  for(i=1;i<=nz;i++){
    cout<<zin[i]<<"\t"<<store[1][i]<<"\t"<<store[2][i]<<"\t"<<store[3][i]<<"\t"<<store[4][i]<<"\t"<<store[5][i]<<"\t"<<store[6][i]<<"\t"<<store[7][i]<<"\t"<<store[8][i]<<endl;
    fout<<zin[i]<<"\t"<<store[1][i]<<"\t"<<store[2][i]<<"\t"<<store[3][i]<<"\t"<<store[4][i]<<"\t"<<store[5][i]<<"\t"<<store[6][i]<<"\t"<<store[7][i]<<"\t"<<store[8][i]<<endl;
  }
  fout.close();

  free_dvector(zin,1,nz);
  free_dmatrix(store,1,8,1,nz);
  */
  ///////////////////////////////////////////////////////////////////////
  // Calculate tauCMB
  //////////////////////////////////////////////////////////////////////
  /*
  file="full_history_popII_st_SB.dat";
  file="full_history.dat";
  int nvar;
  double tau;

  nvar=68;

  tau=a.getTauCMB(file,nvar);

  cout<<tau<<endl;
  */
  ///////////////////////////////////////////////////////////////////////
  // Calculate fraction of mass in collapsed halos 
  //////////////////////////////////////////////////////////////////////
  // Tom Abel was worried that significant fractions of matter would be
  // in collapsed objects that had not formed stars.  This might significantly
  // affect the power spectrum.
  //
  // Defering to him I'll calculate the fraction of mass in collapsed objects
  /*
  double mcool, mjeans, mfilter;
  double fcool, fjeans, ffilter;
  double tjeans;
  double z;

  z=45.0;
  while(z>6.0){
    mcool=coolMass(&c,z);
    mjeans=jeansMass(&c,z);
    mfilter=filterMass(&c,z);
    fcool=c.fColl(z,mcool,0);
    fjeans=c.fColl(z,mjeans,0);
    ffilter=c.fColl(z,mfilter,0);
    tjeans=tvir(&c,mjeans,z,0.6);
    
    cout<<z<<"\t"<<mcool<<"\t"<<mjeans<<"\t"<<mfilter<<"\t"<<fcool<<"\t"<<fjeans<<"\t"<<ffilter<<"\t"<<tjeans<<endl;
    z-=0.5;
  }
  */
  //////////////////////////////////////////////////////////////
  // Test Bubble distribution code - radius
  /////////////////////////////////////////////////////////////
  /*
  double mb,z;
  double Q,R,V,dndlm;
  double rho;

  z=12.0;
  z=zin_arg;

  z=12.0;
  file="bubbledist_z12.dat";
  fout.open(file);
  Q=rz.fillingQ(z);
  cout<<Q<<"\t"<<ZETABUB*c.fColl(z)<<endl;
  cout<<coolMass(&c,z)<<endl;

  fout<<"#"<<Q<<"\t"<<coolMass(&c,z)<<endl;
  rho=CRITDENMSOLMPC*c.getOm0hh();
  R=pow(ZETABUB*coolMass(&c,z)/(4.0*PI*rho/3.0),0.3333);
  while(R<120.0){
    V=4.0*PI*R*R*R/3.0;
    mb=V*rho;
    dndlm=rz.nmBub(mb,z);
    fout<<R<<"\t"<<mb<<"\t"<<dndlm<<"\t"<<3.0*dndlm*V/Q<<endl;
    R*=1.1;
  }
  fout.close();

  z=13.0;
  file="bubbledist_z13.dat";
  fout.open(file);
  Q=rz.fillingQ(z);
  cout<<Q<<"\t"<<ZETABUB*c.fColl(z)<<endl;
  cout<<coolMass(&c,z)<<endl;

  fout<<"#"<<Q<<"\t"<<coolMass(&c,z)<<endl;
  rho=CRITDENMSOLMPC*c.getOm0hh();
  R=pow(ZETABUB*coolMass(&c,z)/(4.0*PI*rho/3.0),0.3333);
  while(R<120.0){
    V=4.0*PI*R*R*R/3.0;
    mb=V*rho;
    dndlm=rz.nmBub(mb,z);
    fout<<R<<"\t"<<mb<<"\t"<<dndlm<<"\t"<<3.0*dndlm*V/Q<<endl;
    R*=1.1;
  }
  fout.close();

  z=14.0;
  file="bubbledist_z14.dat";
  fout.open(file);
  Q=rz.fillingQ(z);
  cout<<Q<<"\t"<<ZETABUB*c.fColl(z)<<endl;
  cout<<coolMass(&c,z)<<endl;

  fout<<"#"<<Q<<"\t"<<coolMass(&c,z)<<endl;
  rho=CRITDENMSOLMPC*c.getOm0hh();
  R=pow(ZETABUB*coolMass(&c,z)/(4.0*PI*rho/3.0),0.3333);
  while(R<120.0){
    V=4.0*PI*R*R*R/3.0;
    mb=V*rho;
    dndlm=rz.nmBub(mb,z);
    fout<<R<<"\t"<<mb<<"\t"<<dndlm<<"\t"<<3.0*dndlm*V/Q<<endl;
    R*=1.1;
  }
  fout.close();

  z=16.0;
  file="bubbledist_z16.dat";
  fout.open(file);
  Q=rz.fillingQ(z);
  cout<<Q<<endl;
  cout<<coolMass(&c,z)<<endl;

  fout<<"#"<<Q<<"\t"<<coolMass(&c,z)<<endl;
  rho=CRITDENMSOLMPC*c.getOm0hh();
  R=pow(ZETABUB*coolMass(&c,z)/(4.0*PI*rho/3.0),0.3333);
  while(R<120.0){
    V=4.0*PI*R*R*R/3.0;
    mb=V*rho;
    dndlm=rz.nmBub(mb,z);
    fout<<R<<"\t"<<mb<<"\t"<<dndlm<<"\t"<<3.0*dndlm*V/Q<<endl;
    R*=1.1;
  }
  fout.close();

  z=18.0;
  file="bubbledist_z18.dat";
  fout.open(file);
  Q=rz.fillingQ(z);
  cout<<Q<<"\t"<<ZETABUB*c.fColl(z)<<endl;
  cout<<coolMass(&c,z)<<endl;

  fout<<"#"<<Q<<"\t"<<coolMass(&c,z)<<endl;
  rho=CRITDENMSOLMPC*c.getOm0hh();
  R=pow(ZETABUB*coolMass(&c,z)/(4.0*PI*rho/3.0),0.3333);
  while(R<120.0){
    V=4.0*PI*R*R*R/3.0;
    mb=V*rho;
    dndlm=rz.nmBub(mb,z);
    fout<<R<<"\t"<<mb<<"\t"<<dndlm<<"\t"<<3.0*dndlm*V/Q<<endl;
    R*=1.1;
  }
  fout.close();
  */
  //////////////////////////////////////////////////////////////
  // average xi code
  /////////////////////////////////////////////////////////////
  /*
  //Integral over rs should give volume when xi_dd=1
  //Integral over bubble mass should then give \bar{b}\bar{Q}
  double r, z, mb;
  double rb,vb,Q;
  z=15.0;
  r=0.01;

  mb=10.0*ZETABUB*coolMass(&c,z);

  vb=mb/(CRITDENMSOLMPC*c.getOm0hh());
  rb=pow(vb/(4.0*PI/3.0),1.0/3.0); 

  Q=rz.fillingQ(z);
  cout<<ZETABUB*c.fColl(z)<<"\t"<<ZETABUB*c.fColl(z,coolMass(&c,z),1)<<"\t"<<ZETABUB*c.fCollPSExact(z,coolMass(&c,z))<<endl;

  cout<<haloBias(mb,z,&c)<<"\t"<<biasm(mb,z,&c)<<"\t"<<biasmST(mb,z,&c)<<"\t"<<c.biasPS(z,mb)<<endl;
  xidd(r,z,&c,1);
  cout<<rz.xiAverageInt(mb,r,z)<<"\t"<<vb<<"\t"<<xidd(2.0,z,&c,0)<<"\t"<<xidd(r,z,&c,0)<<endl;
  //cout<<rz.xiAverageInt(mb,r,z)*rz.meanBiasBub(z)*rz.fillingQ(z)<<endl;
  cout<<rz.meanXiXD(r,z)<<"\t"<<rz.meanBiasBub(z)<<"\t"<<rz.fillingQ(z)<<endl; 
  cout<<rz.meanXiXD(r,z)<<"\t"<<exp(-rz.meanXiXD(r,z))<<endl;
  cout<<-log(1.0-Q)/Q<<"\t"<<Q<<endl;
  cout<<biasmST(mb/400.0,z,&c)<<"\t"<<biasmST(mb/40.0,z,&c)<<"\t"<<biasmST(mb,z,&c)<<endl;

   cout<<rz.correlateXDIntegral(r,z)<<endl;
  cout<<exp(-xidd(rb,z,&c,0)*rz.meanBiasBub(z)*rz.fillingQ(z)*5.0)<<endl;
  */
  //////////////////////////////////////////////////////////////
  // Test xiAverageInt code
  /////////////////////////////////////////////////////////////
  /*
  double r,z,mb,mMin;
  double xi;
  double vb,rb;

  r=1000.0;
  z=7.4;
  rz.setZeta(27.02);
  xidd(r,z,&c,1);

  r=10.0;
  while(r<2500.0){
    cout<<r<<"\t"<<  xidd(r,z,&c,0)<<endl;
    r*=1.1;
  }

  mMin=rz.getZeta()*coolMass(&c,z);
  cout<<mMin<<endl;
  mb=mMin;
  file="xiInt.dat";
  fout.open(file);
  while(mb<1.0e10*mMin){
  vb=mb/(CRITDENMSOLMPC*c.getOm0hh());
  rb=pow(vb/(4.0*PI/3.0),1.0/3.0);
    xi=rz.xiAverageInt(mb,r,z);
    cout<<mb<<"\t"<<rb<<"\t"<<xi<<"\t"<<xi/vb<<endl;
    fout<<mb<<"\t"<<rb<<"\t"<<xi<<"\t"<<xi/vb<<endl;
    mb*=2.0;
  }
  fout.close();

  
  */
  //////////////////////////////////////////////////////////////
  // Test correlation function code
  /////////////////////////////////////////////////////////////
  /*
  double r, z;
  double corrXX,corrDD,corrXD,corrTb;
  z=13.0;

  file="corr.dat";
  fout.open(file);
  fout<<"#"<<z<<"\t"<<rz.fillingQ(z)<<endl;

  cout<<c.fColl(z)<<"\t"<<c.fColl(z,coolMass(&c,z),1)<<endl;

  r=0.01;

 
  //  corrDD=xidd(r,z,&c,1);
  while(r<100.0){
    corrXX=rz.correlateXX(r,z);
    corrXD=rz.correlateXD(r,z);
    corrDD=xidd(r,z,&c,0);
    corrTb=rz.correlateTb(r,z);

    cout<<r<<"\t"<<corrXX<<"\t"<<corrXD<<"\t"<<corrDD<<"\t"<<corrTb<<endl;
    fout<<r<<"\t"<<corrXX<<"\t"<<corrXD<<"\t"<<corrDD<<"\t"<<corrTb<<endl;

    r*=2.0;
  }

  fout.close();
  */
  //////////////////////////////////////////////////////////////
  // Process rawdata file to get simulation x_H history in useful format
  /////////////////////////////////////////////////////////////
  /*
  double *zz, *xH;
  int ndata(91);
  double datum;
  int nz,nx;

  zz=dvector(1,ndata);
  xH=dvector(1,ndata);

  file="rawdata.txt";
  fin.open(file);
  i=1;
  nz=0;
  nx=0;
  while(i<183){
    fin>>datum;
    if(datum>1.0){
      nz++;
      zz[nz]=datum;
    }
    if(datum<1.0){
      nx++;
      xH[nx]=datum;
    }
    i++;
  }
  fin.close();
  
  file="simXH.dat";
  fout.open(file);
  for(i=1;i<=ndata;i++){
    fout<<zz[i]<<"\t"<<xH[i]<<"\t"<<1.0-xH[i]<<endl;
  }
  fout.close();
  */
  //////////////////////////////////////////////////////////////
  // Process rawdata file to get simulation SFR in useful format
  /////////////////////////////////////////////////////////////
  /*
  double *zz, *popII,*popIII;
  int ndata(91);
  double datum;
  int nz,nx;

  zz=dvector(1,ndata);
  popII=dvector(1,ndata);
  popIII=dvector(1,ndata);

  file="sfr_raw.dat";
  fin.open(file);
  for(i=1;i<=ndata;i++){
    fin>>datum;
    zz[i]=datum;
  }
  for(i=1;i<=ndata;i++){
    fin>>datum;
    popII[i]=datum;
  }
  for(i=1;i<=ndata;i++){
    fin>>datum;
    popIII[i]=datum;
  }
  fin.close();
  
  file="simSFR.dat";
  fout.open(file);
  fout<<"#\t z \t popIII \t popII\t [M_sol/Mpc^3/yr]"<<endl;
  for(i=1;i<=ndata;i++){
    fout<<zz[i]<<"\t"<<popII[i]<<"\t"<<popIII[i]<<endl;
  }
  fout.close();
  */
  /////////////////////////////////////////////////////////////
  // Generate ionization history
  ////////////////////////////////////////////////////////////
  /*
  double *zz,*xH,*xI;
  int ndata(91);
  double z, Q;
  double fcoll,fcollST;
  double zeta,zetaST;

  zz=dvector(1,ndata);
  xH=dvector(1,ndata);
  xI=dvector(1,ndata);

  //read in simulation data
  file="simXH.dat";
  fin.open(file);
  for(i=1;i<=ndata;i++){
    fin>>zz[i]>>xH[i]>>xI[i];
  }
  fin.close();

  file="zetaSim.dat";
  fout.open(file);
  fout<<s8<<"\t"<<omb<<"\t"<<om0<<"\t"<<lam0<<"\t"<<h<<"\t"<<n<<endl;
  for(i=1;i<=ndata;i++){
    z=zz[i];
    fcoll=c.fColl(z);
    fcollST=c.fColl(z,coolMass(&c,z),1);
    zeta=xI[i]/fcoll;
    zetaST=xI[i]/fcollST;
    fout<<z<<"\t"<<xH[i]<<"\t"<<zeta<<"\t"<<zetaST<<endl;
    cout<<z<<"\t"<<xH[i]<<"\t"<<zeta<<"\t"<<zetaST<<endl;
  }
  fout.close();
  */
  //////////////////////////////////////////////////////////////
  // Calculate P(k) during reionization
  /////////////////////////////////////////////////////////////
  /*
  double z;
  double corrXX,corrDD,corrXD,corrTb;
  int nr(30), nk(50);
  double *xivec,*pvec,*rvec,*kvec;
  double *xivecXX, *pvecXX;
  double *xivecXD, *pvecXD;
  double *xivecDD, *pvecDD;
  double rstep,kstep;
  double r,k,rmin,rmax,kmin,kmax;
  double ps;
  double xH;

  xivec=dvector(1,nr);
  rvec=dvector(1,nr);
  pvec=dvector(1,nk);
  kvec=dvector(1,nk);

  xivecXX=dvector(1,nr);
  pvecXX=dvector(1,nk);
  xivecXD=dvector(1,nr);
  pvecXD=dvector(1,nk);
  xivecDD=dvector(1,nr);
  pvecDD=dvector(1,nk);

  z=13.0;
  z=zin_arg;
  cout<<z<<endl;
  xidd(r,z,&c,1);

  rmin=0.001;
  rmax=300.0;

  kmin=0.001;
  kmax=100.0;

  //initialise r and k vectors
  rstep=exp(log(rmax/rmin)/(double)(nr-1));
  r=rmin;
  for(i=1;i<=nr;i++){
    rvec[i]=r;
    r*=rstep;
  }
  kstep=exp(log(kmax/kmin)/(double)(nk-1));
  k=kmin;
  for(i=1;i<=nk;i++){
    kvec[i]=k;
    k*=kstep;
  }

  file="corr.dat";
  fout.open(file);
  xH=1.0-rz.fillingQ(z);
  fout<<"#"<<z<<"\t"<<xH<<endl;

  for(i=1;i<=nr;i++){
    r=rvec[i];

    corrXX=rz.correlateXX(r,z);
    corrXD=rz.correlateXD(r,z);
    corrDD=xidd(r,z,&c,0);

    corrTb=corrXX*(1.0+corrDD);
    corrTb+= xH*xH*corrDD;
    corrTb+= corrXD*(2.0*xH+corrXD);

    cout<<r<<"\t"<<corrXX<<"\t"<<corrXD<<"\t"<<corrDD<<"\t"<<corrTb<<endl;
    fout<<r<<"\t"<<corrXX<<"\t"<<corrXD<<"\t"<<corrDD<<"\t"<<corrTb<<endl;

    xivec[i]=corrTb;
    xivecXX[i]=corrXX;
    xivecXD[i]=corrXD;
    xivecDD[i]=corrDD;
  }

  fout.close();

  //Having calculated correlation function take FT to get power spectrum
  cout<<"Taking forrier transform"<<endl;
  FTTableJP(xivec,rvec,nr,pvec,kvec,nk);
  FTTableJP(xivecXX,rvec,nr,pvecXX,kvec,nk);
  FTTableJP(xivecXD,rvec,nr,pvecXD,kvec,nk);
  FTTableJP(xivecDD,rvec,nr,pvecDD,kvec,nk);

  file="ps.dat";
  fout.open(file);

  for(i=1;i<=nk;i++){
    k=kvec[i];
    ps=pvec[i];
    fout<<k<<"\t"<<ps<<"\t"<<pvecXX[i]<<"\t"<<pvecXD[i]<<"\t"<<xH*xH*pvecDD[i]<<"\t"<<pvecDD[i]<<endl;
    cout<<k<<"\t"<<ps<<"\t"<<pvecXX[i]<<"\t"<<pvecXD[i]<<"\t"<<xH*xH*pvecDD[i]<<"\t"<<pvecDD[i]<<endl;
  }
  fout.close();

  free_dvector(rvec,1,nr);
  free_dvector(xivec,1,nr);
  free_dvector(kvec,1,nk);
  free_dvector(pvec,1,nk);

  free_dvector(xivecXX,1,nr);
  free_dvector(pvecXX,1,nk);
  free_dvector(xivecXD,1,nr);
  free_dvector(pvecXD,1,nk);
  free_dvector(xivecDD,1,nr);
  free_dvector(pvecDD,1,nk);
*/
  
  //////////////////////////////////////////////////////////////////////
  // xi average r dependence
  ///////////////////////////////////////////////////////////////////////
  /*
  // Check that <\xi> has correct limiting behaviour... it does
  // Check that <x \delta> has correct limiting behaviour
  //
  double r,z;
  double xi, xi2,xi3;
  double Q,bias;
  double corrXD, corrL;
  double rb,xiAverage;
  double *ans,*result;
  double m, m0;
  double vb;
  z=10.0;
  Q=rz.fillingQ(z);
  bias=rz.meanBiasBub(z);
  rb=rz.meanRBub(z);
  xidd(r,z,&c,1);

  cout<<rb<<"\t"<<bias<<"\t"<<rz.getZeta()<<endl;
  
  ans=dvector(1,1);
  result=dvector(1,1);

  r=100.1;
  setXiXDIntegrand(0.0,ans,ans,z,r,&rz,1);

  m0=rz.getZeta()*coolMass(&c,z);
  m=m0;
  while(m<m0*1e8){
    xiXDIntegrand(m,ans,result);
    vb=m/(CRITDENMSOLMPC*c.getOm0hh());
    rb=pow(vb/(4.0*PI/3.0),1.0/3.0); 
    
    cout<<m<<"\t"<<rb<<"\t"<<result[1]<<"\t"<<rz.nmBub(m,z)<<"\t"<<rz.biasBub(m,z)<<endl;
    m*=2.0;
  }

  
  file="meanXi.dat";
  fout.open(file);

  r=0.001;
  xi=rz.meanXiXD(r,z);
  xi2=rz.meanXiXD(2.0,z);
  while(r<2000.0){
    xi=rz.meanXiXD(r,z);
    xi3=xidd(r,z,&c,0);
    xi2=xidd(rb,z,&c,0);
    xiAverage=rz.meanXiXD(r,z);
    //    xi2=1.0;
    cout<<r<<"\t"<<xi<<"\t"<<Q*bias*xi3<<"\t"<<xidd(r,z,&c,0)<<"\t"<<Q*bias*xi2<<"\t"<<xiAverage<<endl;    
    fout<<r<<"\t"<<xi<<"\t"<<Q*bias*xi3<<"\t"<<xi3<<"\t"<<Q*bias*xi2<<"\t"<<xiAverage<<endl;
    r*=1.5;
  }
  fout.close();
  
  file="meanXiXD.dat";
  fout.open(file);

  r=0.1;
  while(r<1000.0){
    corrXD=rz.correlateXD(r,z);
    xi3=xidd(r,z,&c,0);
    corrL=(1.0-Q)*Q*bias*xi3;  // ignore mean halo bias
    cout<<r<<"\t"<<corrXD<<"\t"<<-corrL<<"\t"<<corrXD/corrL<<endl;    
    fout<<r<<"\t"<<corrXD<<"\t"<<-corrL<<"\t"<<corrXD/corrL<<endl;
    r*=1.5;
  }
  fout.close();
  */
  //////////////////////////////////////////////////////////////////////
  // check bubble overlap code
  ///////////////////////////////////////////////////////////////////////
  /*
  double r,z;
  double xi, xi2,xi3;
  double Q,bias,overlap,overlape;
  double corrXD, corrL;
  double rb;
  z=12.0;
  Q=rz.fillingQ(z);
  bias=rz.meanBiasBub(z);
    rb=rz.meanRBub(z);
  xidd(r,z,&c,1);

  cout<<Q<<"\t"<<rb<<endl;

  file="overlap.dat";
  fout.open(file);

  r=0.1;
  while(r<5000.0){
    overlap=rz.meanBubOverlap(r,z);
    overlape=rz.meanBubOverlapCluster(r,z);
    overlap=  xidd(r,z,&c,0);
    cout<<r<<"\t"<<overlap<<"\t"<<overlape<<"\t"<<(overlap-overlape)/overlap<<"\t"<<Q-overlap<<endl;
    fout<<r<<"\t"<<overlap<<"\t"<<overlape<<"\t"<<(overlap-overlape)/overlap<<"\t"<<Q-overlap<<endl;
    r*=1.3;
  }
  fout.close();
 */
  ////////////////////////////////////////////////////////////
  // Renormalisation factor
  ////////////////////////////////////////////////////////////
  /*
  double Q,renorm;
  Q=0.0;
  while(Q<1.1){
    Q+=0.01;
    renorm= -log(1.0-Q)/Q;
    cout<<Q<<"\t"<<renorm<<endl;
  }
*/
  ////////////////////////////////////////////////////////////
  // Plot XD integral kernel
  ////////////////////////////////////////////////////////////
  /*
  double z,mmin,m;
  double Km,r;
  double *ans;

  ans=dvector(1,1);

  z=15.0;
  r=0.001;
  xidd(r,z,&c,1);
  rz.setZeta(40.0);
  setXiXDIntegrand(0.0,ans,ans,z,r,&rz,1);

  mmin=rz.getZeta()*coolMass(&c,z);

  m=mmin;

  file="XDkernel.dat";
  fout.open(file);
  while(m<10000.0*mmin){
    
    xiXDIntegrand(m,ans,ans);
    Km=ans[1];
    cout<<m<<"\t"<<Km<<endl;
    fout<<m<<"\t"<<Km<<endl;
    m*=2.0;
  }
  fout.close();

  free_dvector(ans,1,1);
  */
  //////////////////////////////////////////////////////////////
  // Calculate bubble distribution for simulation
  /////////////////////////////////////////////////////////////
  /*
  double mb,z;
  double Q,R,V,dndlm;
  double rho;
  int nz(7);
  double *zlist;

  double *zzSim,*zetaSim, trash;
  int ndataSim(91);
  char tag[20];
  Spline zetaSP;


  zlist=dvector(1,nz);
  zlist[1]=24.8;
  zlist[2]=16.6;
  zlist[3]=12.8;
  zlist[4]=10.6;
  zlist[5]=8.1;
  zlist[6]=6.6;
  zlist[7]=6.0;


  //calculate zeta value for the specified redshift
  zzSim=dvector(1,ndataSim);
  zetaSim=dvector(1,ndataSim);

  file="zetaSim.dat";
  fin.open(file);
  for(i=1;i<=ndataSim;i++){
    fin>>zzSim[i]>>trash>>zetaSim[i]>>trash;
    //  cout<<zzSim[i]<<"\t"<<zetaSim[i]<<endl;
  }
  fin.close();

  zetaSP.setSplineSP(ndataSim,zzSim,zetaSim);

  //Calculate bubble distribution at each redshift

  for(i=1;i<=nz;i++){
    z=zlist[i];
    sprintf(tag,"bubbledist_z_%05.2lf.dat",z);
    fout.open(tag);
    cout<<z<<endl;
    rz.setZeta(zetaSP.returnValue(z));
    cout<<rz.getZeta()<<endl;
    Q=rz.fillingQ(z);
    cout<<z<<"\t"<<Q<<"\t"<<rz.getZeta()*c.fColl(z)<<"\t"<<rz.getZeta()<<endl;
    cout<<coolMass(&c,z)<<endl;
    
    fout<<"#"<<Q<<"\t"<<coolMass(&c,z)<<"\t"<<rz.getZeta()<<endl;
    rho=CRITDENMSOLMPC*c.getOm0hh();
    R=pow(rz.getZeta()*coolMass(&c,z)/(4.0*PI*rho/3.0),0.3333);
    while(R<120.0){
      V=4.0*PI*R*R*R/3.0;
      mb=V*rho;
      dndlm=rz.nmBub(mb,z);
      //  cout<<R<<"\t"<<mb<<"\t"<<dndlm<<"\t"<<3.0*dndlm*V/Q<<endl;
      fout<<R<<"\t"<<mb<<"\t"<<dndlm<<"\t"<<3.0*dndlm*V/Q<<endl;
      R*=1.1;
    }
    fout.close();

  }

  free_dvector(zzSim,1,ndataSim);
  free_dvector(zetaSim,1,ndataSim);
  free_dvector(zlist,1,nz);
  */
  //////////////////////////////////////////////////////////////
  // Calculate P(k) during reionization - using simulation zeta
  /////////////////////////////////////////////////////////////
  /*
  double z;
  double corrXX,corrDD,corrXD,corrTb;
  int nr, nk(50), shift(30),shift2(35),shift3(20);
  double *xivec,*pvec,*rvec,*kvec;
  double *xivecXX, *pvecXX;
  double *xivecXD, *pvecXD;
  double *xivecDD, *pvecDD;
  double rstep,kstep;
  double r,k,rmin,rmax,kmin,kmax;
  double ps;
  double xH;
  int fflag;
  double *zzSim,*zetaSim, trash;
  int ndataSim(91);
  double tbright,tk,lya,xi;
  double Q,renorm;
  double s8t,ombt,om0t,lam0t,ht,nt;
  char tag[50];
  Spline zetaSP;

  //    shift=90;
  //shift2=180;
  //shift3=90;
    shift=20;
   shift2=40;
   shift3=20;

  nr=shift+shift2+shift3;
  nr--;
  nk=150;

  //  nr=300;
  nk=500;

  xivec=dvector(1,nr);
  rvec=dvector(1,nr);
  pvec=dvector(1,nk);
  kvec=dvector(1,nk);

  xivecXX=dvector(1,nr);
  pvecXX=dvector(1,nk);
  xivecXD=dvector(1,nr);
  pvecXD=dvector(1,nk);
  xivecDD=dvector(1,nr);
  pvecDD=dvector(1,nk);

  z=zin_arg;
  cout<<z<<endl;
  xidd(r,z,&c,1);

  //calculate zeta value for the specified redshift
  zzSim=dvector(1,ndataSim);
  zetaSim=dvector(1,ndataSim);

  file="zetaSim.dat";
  fin.open(file);
  fin>>s8t>>ombt>>om0t>>lam0t>>ht>>nt;

  if( fabs(s8-s8t)>1.0e-4 ||  fabs(omb-ombt)>1.0e-4  ||  fabs(om0-om0t)>1.0e-4
       ||  fabs(lam0-lam0t)>1.0e-4 ||  fabs(h-ht)>1.0e-4 
      ||  fabs(n-nt)>1.0e-4){
    cout<<"cosmological parameter mismatch for zeta spline"<<endl;
  }

  for(i=1;i<=ndataSim;i++){
    fin>>zzSim[i]>>trash>>zetaSim[i]>>trash;
  }
  fin.close();

  zetaSP.setSplineSP(ndataSim,zzSim,zetaSim);
  rz.setZeta(zetaSP.returnValue(z));

  free_dvector(zzSim,1,ndataSim);
  free_dvector(zetaSim,1,ndataSim);

  //  rz.setZeta(40.0);
  cout<<rz.getZeta()<<endl;


  xH=1.0-rz.fillingQ(z);
  cout<<z<<"\t"<<rz.getZeta()<<"\t"<<xH<<"\t"<<rz.meanBiasBub(z)<<"\t"<<rz.meanRBub(z)<<endl;

  cout<<"xH="<<xH<<"\t mean Rbub="<<rz.meanRBub(z)<<endl;

  //calculate mean brightness temperature at this redshift
  xi=1.0-xH;
  lya=0.0;       //need prescription for lya and tk
  tk=0.0;
  tbright=tocm.tBrightGen(z,tk,xi,lya);

  //limits of integration
  
  rmin=1.0e-4;
  //  rmin=1.0e-3;
  //rmin=1.0e-5;

  kmin=0.001;
  kmax=100.0;
  //  kmax=10.0;
  //  kmin=0.01;
  //  kmax=30.0;
   //   kmax=10.0;
   kmax=30.0;
  //kmax=20.0;
  //  kmax=15.0;
  // kmax=10.0;

  //calculate r values

  //  rmax=60.0;
  rmax=10.0;
  //  rmax=1.0;
  //initialise r and k vectors
  rstep=exp(log(rmax/rmin)/(double)(shift-1));
  r=rmin;
  for(i=1;i<=shift;i++){
    rvec[i]=r;
    r*=rstep;
  }
  
  rmin=r;
  rmax=300.0;
  rstep=exp(log(rmax/rmin)/(double)(shift2-1));
  //  rstep=(rmax-r)/(double)(shift-1);
  for(i=shift+1;i<=shift+shift2;i++){
    rvec[i]=r;
    r*=rstep;
    //r+=1.0;
  }

  rmin=r;
  //  rmax=2000.0;
  rmax=2200.0;
  // rmax=2450.0;
  //  rmax=8000.0;
  //  rmax=9500.0;
  //rmax=9200.0;
  //rmax=20000.0;
  rstep=exp(log(rmax/rmin)/(double)(shift3-1));
  //  rstep=(rmax-r)/(double)(shift-1);
  for(i=shift+shift2+1;i<=nr;i++){
    rvec[i]=r;
    r*=rstep;
  }
  
  kmin=0.001;
  kmax=31.0;
  kstep=exp(log(kmax/kmin)/(double)(nk-1));
  k=kmin;
  for(i=1;i<=nk;i++){
    kvec[i]=k;
    k*=kstep;
  }
  

   cout<<z<<"\t"<<rz.getZeta()<<"\t"<<xH<<"\t"<<rz.meanBiasBub(z)<<"\t"<<rz.meanRBub(z)<<endl;
   Q=rz.fillingQ(z);
  renorm= -log(1.0-Q)/Q;
  cout<<"z="<<z<<"\t Q="<<Q<<"\t renorm="<<renorm<<"\t"<<rz.getZeta()*c.fColl(z)<<endl;
   cout<<"bub bias:"<<rz.meanBiasBub(z)<<"\t halo bias:"<<meanHaloBias(z,&c)<<"\t bound voilation:"<<1.0-xH*pow(meanHaloBias(z,&c),2.0)<<endl;  

   fflag=1;
   if(fflag==0){
     cout<<"calculating correlation functions"<<endl;

  sprintf(tag,"corr_z_%05.2lf.dat",z);
  //file="corr.dat";
  fout.open(tag);
  fout.precision(12);
  
  fout<<"#"<<z<<"\t"<<xH<<"\t"<<rz.getZeta()<<endl;

  for(i=1;i<=nr;i++){
    r=rvec[i];

    cout<<endl;
    cout<<"r="<<r<<endl;
 
    cout<<"correlateXX...";
    corrXX=rz.correlateXX(r,z);
    cout<<"correlateXD...";
    corrXD=rz.correlateXD(r,z);
    cout<<"correlateDD...";
    corrDD=xidd(r,z,&c,0);
    cout<<"done"<<endl;

    corrTb=corrXX*(1.0+corrDD);
    corrTb+= xH*xH*corrDD;
    corrTb+= corrXD*(2.0*xH+corrXD);

    cout<<sqrt(fabs(corrXX*corrDD))<<"\t"<<sqrt(corrXD*corrXD)<<"\t"<<sqrt(fabs(corrXX*corrDD))-fabs(corrXD)<<endl;

    cout<<r<<"\t"<<corrXX<<"\t"<<corrXD<<"\t"<<corrDD<<"\t"<<corrTb<<endl;
    fout<<r<<"\t"<<corrXX<<"\t"<<corrXD<<"\t"<<corrDD<<"\t"<<corrTb<<endl;

    xivec[i]=corrTb;
    xivecXX[i]=corrXX;
    xivecXD[i]=-corrXD;
    xivecDD[i]=corrDD;
  }

  fout.close();
  
   }else{
     //read in data from file
     sprintf(tag,"corr_z_%05.2lf.dat",z);
     fin.open(tag);
     fin.getline(tag,50);
     for(i=1;i<=nr;i++){
       fin>>rvec[i]>>xivecXX[i]>>xivecXD[i]>>xivecDD[i]>>xivec[i];
     }
     fin.close();
   }
  
  //Having calculated correlation function take FT to get power spectrum
  cout<<"Taking forrier transform: Tb"<<endl;
  FTTable(xivec,rvec,nr,pvec,kvec,nk);
  cout<<"Taking forrier transform: XX"<<endl;
  FTTable(xivecXX,rvec,nr,pvecXX,kvec,nk);
  cout<<"Taking forrier transform: XD"<<endl;
  FTTable(xivecXD,rvec,nr,pvecXD,kvec,nk); // numerical accuracy poor here
  cout<<"Taking forrier transform: DD"<<endl;
  FTTable(xivecDD,rvec,nr,pvecDD,kvec,nk);


  sprintf(tag,"ps_z_%05.2lf.dat",z);
  //file="ps.dat";
  fout.open(tag);

  for(i=1;i<=nk;i++){
    k=kvec[i];
    ps=pvec[i];
    fout<<k<<"\t"<<ps<<"\t"<<pvecXX[i]<<"\t"<<pvecXD[i]<<"\t"<<xH*xH*pvecDD[i]<<"\t"<<pvecDD[i]<<"\t"<<ps/xH/xH<<endl;
    cout<<k<<"\t"<<ps<<"\t"<<pvecXX[i]<<"\t"<<pvecXD[i]<<"\t"<<xH*xH*pvecDD[i]<<"\t"<<pvecDD[i]<<"\t"<<ps/xH/xH<<endl;
  }
  fout.close();

  free_dvector(rvec,1,nr);
  free_dvector(xivec,1,nr);
  free_dvector(kvec,1,nk);
  free_dvector(pvec,1,nk);

  free_dvector(xivecXX,1,nr);
  free_dvector(pvecXX,1,nk);
  free_dvector(xivecXD,1,nr);
  free_dvector(pvecXD,1,nk);
  free_dvector(xivecDD,1,nr);
  free_dvector(pvecDD,1,nk);
  */

  ////////////////////////////////////////////////////////////////////////
  // Calculate power spectrum D + X + T + A
  ////////////////////////////////////////////////////////////////////////
 
  double *zN, *kN;
  double **betaN, **waN, **gtN, **thermalN;
  double *beta, *thermalI;
  double **powerDD, **powerXD, **powerXX, **powerTB, **powerFZH;
  double **powerPostTB;
  
  double **gtNU;
  double *gtV, **gtM;
  double *rFN, *rTN;
  double *powerDDV, *powerXDV, *powerXXV, *powerTBV;
  double *corrDDV, *corrXDV, *corrXXV, *corrTBV;

  int fast_flag(0);
  //  int nz(120), nk(5);
  //int nz(149), nk(3);
  int nz(150), nk(27);
  //  int nz(20), nk(27);


  int nr(79);
  double *rN;

  double z,k,tk,xi,xe,lya,xH;
  double *tsN;
  double zstart, zend, zstep(0.0),kstart, kend, kstep(0.0);
  char tag[50];

  zN=dvector(1,nz);
  kN=dvector(1,nk);

  betaN=dmatrix(1,nz,1,5);
  waN=dmatrix(1,nz,1,nk);
  gtN=dmatrix(1,nz,1,nk);
  thermalN=dmatrix(1,nz,1,4);
  gtNU=dmatrix(1,nz,1,nk);

  tsN=dvector(1,nz);
  rTN=dvector(1,nz);
  rFN=dvector(1,nz);  
  beta=dvector(1,5);
  thermalI=dvector(1,3);
  gtM=dmatrix(1,4,1,nz);
  gtV=dvector(1,4);

  powerDD=dmatrix(1,nz,1,nk);
  powerXD=dmatrix(1,nz,1,nk);
  powerXX=dmatrix(1,nz,1,nk);
  powerTB=dmatrix(1,nz,1,nk);
  powerFZH=dmatrix(1,nz,1,nk);   //FZH D+X Tb power spectrum
  powerPostTB=dmatrix(1,nz,1,nk);  //post reionization TB power spectrum

  corrDDV=dvector(1,nr);
  corrXDV=dvector(1,nr);
  corrXXV=dvector(1,nr);
  corrTBV=dvector(1,nr);

  powerDDV=dvector(1,nk);
  powerXDV=dvector(1,nk);
  powerXXV=dvector(1,nk);
  powerTBV=dvector(1,nk);

  rN=dvector(1,nr);

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

  a.globalHistory();
  cout<<"tau="<<a.getTauCMB()<<endl;
  //specify z and k values
  //  zstart=8.0;
  //  zstart=3.5;

  //zstart=7.0;
  //  zend=25.0;
  //zend=100.0;
  //  zend=200.0;
  
  zstart=3.5;
  zend=400.0;

  //  zstart=0.1;
  //zend=25.0;

  //  zstart=3.5;
  //zend=500.0;

  //  kstart=0.1;
  //kend=1.0;
  //kstart=0.1;
  //kend=10.0;
  //  kstart=0.03163;
  //kend=3.163;

  kstart=0.001;
  kend=3163.0;

  //  kstart=0.01;
  //kend=100.0;

  // if(nz>1) zstep=(zend-zstart)/(double)(nz-1);
  //  for(i=1;i<=nz;i++){
  //  zN[i]=zstart+(double)(i-1)*zstep;
  //}

  if(nz>1) zstep=log(zend/zstart)/(double)(nz-1);
  for(i=1;i<=nz;i++){
    zN[i]=zstart*exp((double)(i-1)*zstep);
    //    cout<<zN[i]<<endl;
  }

 if(nk>1) kstep=exp(log(kend/kstart)/(double)(nk-1));
  k=kstart;
  for(i=1;i<=nk;i++){
    kN[i]=k;
    //  cout<<kN[i]<<endl;
    k*=kstep;
  }

  // obtain thermal history and beta values at different z
  cout<<"Calculting thermal history"<<endl;
  for(i=1;i<=nz;i++){
    z=zN[i];
    a.getTIGM(z,thermalI);
    tk=thermalI[1];
    xi=thermalI[2];
    xe=thermalI[3];
    lya=a.lyaFlux(z);
    tocm.getBetaTb(z,tk,xe,lya,beta);

    tsN[i]=tocm.tSpinGen(z,tk,xe,lya);
    //store beta values in redshift bin
    for(j=1;j<=5;j++){
      betaN[i][j]=beta[j];
    }
    //store thermal info in redshift bin
    for(j=1;j<=3;j++){
      thermalN[i][j]=thermalI[j];
    }
    thermalN[i][4]=lya;  //save lya data
    //    cout<<zN[i]<<"\t"<<thermalN[i][4]<<endl;
  }

  //writing thermal history
  file="tk_history_z.dat";
  fout.open(file);
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t"<<thermalN[i][1]<<"\t"<<thermalN[i][2]<<"\t"<<thermalN[i][3]<<"\t"<<tsN[i]<<"\t"<<TCMB*(1.0+zN[i])<<"\t"<<thermalN[i][4]<<endl;
  }
  fout.close();

  file="beta_history_z.dat";
  fout.open(file);
  for(i=1;i<=nz;i++){
    fout<<zN[i];
    for(j=1;j<=5;j++){
      fout<<"\t"<<betaN[i][j];
    }
    fout<<endl;
  }
  fout.close();

  //calculating cutoff scales
  cout<<"calculating cutoff scales"<<endl;
  for(i=1;i<=nz;i++){
    rFN[i]=tocm.getCutoffFilter(zN[i]);
    rTN[i]=tocm.getCutoffThermal(zN[i],thermalN[i][1]);
    //    cout<<zN[i]<<"\t"<<rFN[i]<<"\t"<<rTN[i]<<endl;
  }

  file="cutoff_scales_z.dat";
  fout.open(file);
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t"<<rFN[i]<<"\t"<<rTN[i]<<endl;
  }
  fout.close();

  //calculate Lya window function in z and k
  cout<<"Calculating L window function"<<endl;
  for(i=1;i<=nz;i++){
    for(j=1;j<=nk;j++){
      z=zN[i];
      k=kN[j];
      if(fast_flag==0)  waN[i][j]=tocm.getWindowLya(z,k);
    }
  }

  //writing L window function
  file="L_window.dat";
  fout.open(file);
  for(j=1;j<=nk;j++){
    fout<<kN[j];
    for(i=1;i<=nz;i++){
      fout<<"\t"<<waN[i][j];
    }
    fout<<endl;
  }
  fout.close();

  //calculate Temp window function in z and k
  cout<<"Calculating T window function"<<endl;
  if(nz>2){
    for(j=1;j<=nk;j++){
      k=kN[j];
      if(fast_flag==0) tocm.getGTN(zN,k,gtM,nz);
      for(i=1;i<=nz;i++){
	gtN[i][j]=gtM[1][i];
	gtNU[i][j]=gtM[2][i];
      }
    }
  } else{
    for(j=1;j<=nk;j++){
      k=kN[j];
      for(i=1;i<=nz;i++){
	z=zN[i];
	if(fast_flag==0) tocm.getGT(z,k,gtV);  //1:inhomo  2: uniform  3: adiabatic
	gtN[i][j]=gtV[1];
	gtNU[i][j]=gtV[2];
      }
    }
  }

  //writing L window function
  file="T_window.dat";
  fout.open(file);
  for(j=1;j<=nk;j++){
    fout<<kN[j];
    for(i=1;i<=nz;i++){
      fout<<"\t"<<gtN[i][j];
    }
    fout<<endl;
  }
  fout.close();

  //writing L window function
  file="T_window_uniform.dat";
  fout.open(file);
  for(j=1;j<=nk;j++){
    fout<<kN[j];
    for(i=1;i<=nz;i++){
      fout<<"\t"<<gtNU[i][j];
    }
    fout<<endl;
  }
  fout.close();

  //Now calculate Density and Ionization power spectra
  cout<<"calculating D and X power spectra"<<endl;

  //get radial values at which to calculate correlation function
  rz.getRVec(rN);

  for(i=1;i<=nr;i++){
    cout<<rN[i]<<endl;
  }

  cout<<"getting correlation functions"<<endl;
  file="power_store.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(j=1;j<=nk;j++){
    fout<<kN[j]<<"\t";
  }
  fout<<endl;

  for(i=1;i<=nz;i++){
    z=zN[i];
    fout<<z<<"\t";
    if(z>ZSFR || thermalN[i][2]<rz.getXICUT()){
      //no star formation and well into linear regime
      //or very little star formation so ionization can be ignored
      for(j=1;j<=nk;j++){
	powerDDV[j]=Plinear(kN[j],zN[i],&c);
	powerXDV[j]=0.0;
	powerXXV[j]=0.0;
	powerTBV[j]=powerDDV[j];
      }
    }else if(z<a.getZReion()){
      //in fully ionized regime need post-reionization code
      cout<<"DD...";
      rz.getCorrDD(z,rN,corrDDV,nr);
      FTTable(corrDDV,rN,nr,powerDDV,kN,nk);      
      for(j=1;j<=nk;j++){
	powerXDV[j]=0.0;
	powerXXV[j]=0.0;
	powerTBV[j]=0.0;
      }
      cout<<"Done"<<endl;
    }else{   
      //in star forming regime need ionization and non-linear density
      xH=1.0-thermalN[i][2];  //need bubble xH here
      cout<<"DD...";
      rz.getCorrDD(z,rN,corrDDV,nr);
      FTTable(corrDDV,rN,nr,powerDDV,kN,nk);
      cout<<"XD...";
      rz.getCorrXD(z,rN,corrXDV,nr);
      FTTable(corrXDV,rN,nr,powerXDV,kN,nk);
      cout<<"XX...";
      rz.getCorrXX(z,rN,corrXXV,nr);
      FTTable(corrXXV,rN,nr,powerXXV,kN,nk);
      cout<<"TB...";
      rz.getCorrTB(z,rN,corrTBV,nr,corrDDV,corrXDV,corrXXV,xH);
      FTTable(corrTBV,rN,nr,powerTBV,kN,nk);
      cout<<"done"<<endl;
    }

    //handle post reionization power spectrum
    for(j=1;j<=nk;j++){
      powerPostTB[i][j]=tocm.postReionPowerSpectrum(zN[i],kN[j],powerDDV[j]);
    }

    //place vector info into the full matrix
    for(j=1;j<=nk;j++){
      powerDD[i][j]=powerDDV[j];
      powerXD[i][j]=powerXDV[j];
      powerXX[i][j]=powerXXV[j];
      powerFZH[i][j]=powerTBV[j];
      ///////////////////
      // Calculating full power spectrum on the fly
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerFZH[i][j],powerDD[i][j],powerXD[i][j],betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtN[i][j],waN[i][j],thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);

      fout<<powerTB[i][j]<<"\t";
    
    }
    fout<<endl;
  }
  fout.close();


  /////////////////////////////////////////////////////////////////////
  // Output all data for post processing
  // first output all data into multiple files
  for(i=1;i<=nz;i++){
    sprintf(tag,"data/powerTB_z_%05.2lf.dat",zN[i]);
    fout.open(tag);
    fout.precision(12); 
    for(j=1;j<=nk;j++){
      fout<<zN[i]<<"\t"<<kN[j]<<"\t"<<powerTB[i][j]<<"\t"<<powerFZH[i][j]<<"\t"<<powerDD[i][j]<<"\t"<<powerXD[i][j]<<"\t"<<powerXX[i][j]<<"\t"<<betaN[i][1]<<"\t"<<betaN[i][2]<<"\t"<<betaN[i][3]<<"\t"<<betaN[i][4]<<"\t"<<betaN[i][5]<<"\t"<<gtN[i][j]<<"\t"<<gtNU[i][j]<<"\t"<<waN[i][j]<<"\t"<<thermalN[i][1]<<"\t"<<thermalN[i][2]<<"\t"<<thermalN[i][3]<<endl;
    }
    fout.close();
  }

  ////////////////////////////////
  //From above information calculate Tb power spectrum
  cout<<"calculating TB power spectrum"<<endl;

  ///////////////////////////////////////////////////////////////
  //Output power spectrum information for use in plotting
  cout<<"writing Tb data to file"<<endl;
  file="powerTb.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerFZH[i][j],powerDD[i][j],powerXD[i][j],betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtN[i][j],waN[i][j],thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  cout<<"writing FZH data to file"<<endl;

  file="powerTb_FZH.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerFZH[i][j],powerDD[i][j],powerXD[i][j],betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],0.0,0.0,thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<endl;

      //      fout<<powerFZH[i][j]*betaN[i][1]*betaN[i][1]+(2.0/3.0*betaN[i][1]*betaN[i][5]+1.0/5.0*betaN[i][5]*betaN[i][5])*powerDD[i][j]+2.0/3.0*betaN[i][2]*betaN[i][5]*powerXD[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  file="powerTb_DD.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerDD[i][j],powerDD[i][j],powerXD[i][j],betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],0.0,0.0,thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<endl;
      //      fout<<powerDD[i][j]*(betaN[i][1]*betaN[i][1]+2.0/3.0*betaN[i][1]*betaN[i][5]+1.0/5.0*betaN[i][5]*betaN[i][5])<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  cout<<"writing FZH data to file"<<endl;

  file="powerFZH.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      fout<<powerFZH[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  cout<<"writing DD data to file"<<endl;

  file="powerDD.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      fout<<powerDD[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  cout<<"writing XX data to file"<<endl;

  file="powerXX.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      fout<<powerXX[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  cout<<"writing XD data to file"<<endl;

  file="powerXD.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      fout<<powerXD[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  cout<<"writing Post Reionization data to file"<<endl;

  file="powerPostTB.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      fout<<powerPostTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();


  //special part : uniform heating

  cout<<"calculating TB power spectrum: uniform heating"<<endl;
  //Output power spectrum information
  cout<<"writing Tb data to file: uniform heating"<<endl;
  file="powerTb_uniform.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerFZH[i][j],powerDD[i][j],powerXD[i][j],betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtNU[i][j],waN[i][j],thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  //special part: no ionization
  cout<<"calculating TB power spectrum: no ionization"<<endl;
  //Output power spectrum information
  cout<<"writing Tb data to file: no ionization"<<endl;
  file="powerTb_noion.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerDD[i][j],powerDD[i][j],0.0,betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtN[i][j],waN[i][j],thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  //no lya

  cout<<"calculating TB power spectrum: no lya"<<endl;

  //Output power spectrum information
  cout<<"writing Tb data to file: no lya"<<endl;
  file="powerTb_nolya.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerFZH[i][j],powerDD[i][j],powerXD[i][j],betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtN[i][j],0.0,thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  //ionization only
  cout<<"calculating TB power spectrum: ionization only"<<endl;
  //Output power spectrum information
  cout<<"writing Tb data to file: ionization only"<<endl;
  file="powerTb_ion.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerFZH[i][j],powerDD[i][j],powerXD[i][j],betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtNU[i][j],0.0,thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();


  //lya only
  cout<<"calculating TB power spectrum: lya only"<<endl;
  //Output power spectrum information
  cout<<"writing Tb data to file: lya only"<<endl;
  file="powerTb_lya.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerDD[i][j],powerDD[i][j],0.0,betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtNU[i][j],waN[i][j],thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();


  //temp only
  cout<<"calculating TB power spectrum: temp only"<<endl;
  //Output power spectrum information
  cout<<"writing Tb data to file: temp only"<<endl;
  file="powerTb_temp.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerDD[i][j],powerDD[i][j],0.0,betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtN[i][j],0.0,thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  //desnity only
  cout<<"calculating TB power spectrum: density only"<<endl;
  //Output power spectrum information
  cout<<"writing Tb data to file: density only"<<endl;
  file="powerTb_dens.dat";
  fout.open(file);
  fout<<0.0<<"\t";
  for(i=1;i<=nk;i++){
    fout<<kN[i]<<"\t";
  }
  fout<<endl;
  for(i=1;i<=nz;i++){
    fout<<zN[i]<<"\t";
    for(j=1;j<=nk;j++){
      powerTB[i][j]=tocm.fullPowerSpectrum(zN[i],kN[j],powerDD[i][j],powerDD[i][j],0.0,betaN[i][1],betaN[i][2],betaN[i][3],betaN[i][4],betaN[i][5],gtNU[i][j],0.0,thermalN[i][1], thermalN[i][2], rFN[i], rTN[i], powerPostTB[i][j]);
      fout<<powerTB[i][j]<<"\t";
    }
    fout<<endl;
  }
  fout.close();

  ////////////////////////////////
  // Save parameter set to file



  free_dvector(tsN,1,nz);

  ///////////////////////////////////////////////////////////////////////
  return 0;
}




