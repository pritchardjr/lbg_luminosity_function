
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
#include "galaxies.h"
#include "fisherLF.h"
#include "astronomy.h"

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
  string files;
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

/////////////////////////////////
// dcosmology testing
/////////////////////////////////
  cout<<"testing"<<endl;
  cout<<c.getScale()<<endl;
  cout<<c.powerSpectrum(0.1)<<endl; 
  double dummy, k;
  //cout<<c.TFMaster(0.1,dummy,0)<<endl;
  c.resetPowerSpectrum(); 
  cout<<setSigmatop(0.1,&c,1)<<endl;
  cout<<sigmatop(0.1)<<endl;
  cout<<c.sigma0fM(1.0e8,dummy,0)<<endl;
  cout<<c.biasPS(1.0,1.0e8)<<endl;
  cout<<c.dndlM(1.0,1.0e8)<<endl;
  k=0.001;
  while(k<10.0){
     cout<<k<<"\t"<<c.TFMaster(k,dummy,0)<<endl;
     k*=10.0;
  }

////////////////////////////////////////////////////////////////////////
// Test code for getting mfp from gamma
///////////////////////////////////////////////////////////////////////
/*
double mfp, gamma, z;
string files, tag;

Fisher ff;

z=1.0;
for(i=2;i<=6;i++){
	z+=1.0;
	cout<<z<<endl;
	a.initPVParamMEHR(z);

     tag="mfp_";
     files=tag+ff.numToString(z)+".dat";
     cout<<files<<endl;
     fout.open(files.c_str());

	gamma=0.01;
	while(gamma<100.0){
	cout<<gamma<<"\t"<<lyf.deltaIonFromGamma(z,gamma)<<"\t"<<lyf.mfpFromGamma(z,gamma)<<endl;
	fout<<gamma<<"\t"<<lyf.deltaIonFromGamma(z,gamma)<<"\t"<<lyf.mfpFromGamma(z,gamma)<<endl;

	gamma*=1.2;
	}

	fout.close();
}
*/
////////////////////////////////////////////////////////////////////////
// Test code for calculating ionization history
///////////////////////////////////////////////////////////////////////
/*
string files;
double Nion, xi, tau,z, gamma, zeta, mfp, clump;
vector<double> param;
double zstep(0.15);

files="output.dat";
fout.open(files.c_str());

Ionization Ion(&c,&a);
param.clear();
param.push_back(1);
param.push_back(45.0);
param.push_back(0.0);
param.push_back(0.0);
param.push_back(0.0);

//two step param
//param.clear();
//param.push_back(2);
//param.push_back(45.0);
//param.push_back(120.0);
//param.push_back(13.0);
//param.push_back(2.0);

//Nion param
//param.clear();
//param.push_back(3);
//param.push_back(1.0);
//param.push_back(15.0);
//param.push_back(0.0);
//param.push_back(0.0);
//param.push_back(0.0);


Ion.setZetaParam(0.0,param,1);
cout<<Ion.getTauCMB()<<endl;
fout<<Ion.getTauCMB()<<endl;

z=20.0;
zstep=0.5;
while(z>2.5){
	Nion=Ion.getNion(z)*MPC*MPC*MPC;
	xi=Ion.getXI(z);
	gamma=lyf.getGammaFromNdot(z,Nion);
	zeta=Ion.getZeta(z);
	mfp=lyf.mfpFromGamma(z,gamma);
	a.initPVParamMEHR(z);
	clump=rz.clumpingBub(z,xi);
//	clump=4.0*a.clumpingMEHR(z,xi);
	cout<<z<<"\t"<<Nion<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<Ion.getClumping(z,xi)<<"\t"<<mfp<<"\t"<<clump<<endl;
	fout<<z<<"\t"<<Nion<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<Ion.getClumping(z,xi)<<"\t"<<mfp<<"\t"<<clump<<endl;	
	z-=zstep;
}

fout.close();

*/
////////////////////////////////////////////////////////////////////////
// Grid likelihood values
///////////////////////////////////////////////////////////////////////
/*
string files;
double zeta0, zeta1;
double like;
vector<double> param;
Ionization Ion(&c,&a);

files="likelihoodZeta.dat";
fout.open(files.c_str());

for(i=1;i<=21;i++){
	zeta0=30+(double)(i-1)*2.5;
	for(j=1;j<=21;j++){
		zeta1=0.0+(double)(j-4)*0.01;

		param.clear();
		param.push_back(zeta0);
		param.push_back(zeta1);
		Ion.setZetaParam(0.0,param,1);
		like=Ion.likelihoodZeta(param);
		cout<<zeta0<<"\t"<<zeta1<<"\t"<<like<<"\t"<<Ion.getXI(6.0)<<"\t"<<Ion.getXI(8.0)<<"\t"<<Ion.getXI(10.0)<<endl;
		fout<<zeta0<<"\t"<<zeta1<<"\t"<<like<<"\t"<<Ion.getXI(6.0)<<"\t"<<Ion.getXI(8.0)<<"\t"<<Ion.getXI(10.0)<<endl;
	}
}
fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Grid likelihood values - 1 parameter constant zeta
///////////////////////////////////////////////////////////////////////
/*
string files;
double zeta0,z;
double like;
vector<double> param;
Ionization Ion(&c,&a);
int k;
int n1(201);
double p1min,p1max,p1step;

files="likelihoodZeta_const_zeta.dat";
fout.open(files.c_str());

p1min=30.0;
p1max=70.0;

//n1=301;
//p1min=10.0;
//p1max=300.0;

p1step=(p1max-p1min)/(double)(n1-1);

for(i=1;i<=n1;i++){
	zeta0=p1min+p1step*(double)(i-1);
		
	param.clear();
	param.push_back(1);
	param.push_back(zeta0);
	param.push_back(0.0);
	param.push_back(0.0);
	param.push_back(0.0);

		Ion.setZetaParam(0.0,param,1);
		like=Ion.likelihoodZeta(param);
		cout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<like<<"\t"<<Ion.getXI(6.0)<<"\t"<<Ion.getXI(7.0)<<"\t"<<Ion.getXI(8.0)<<"\t"<<Ion.getXI(9.0)<<"\t"<<Ion.getXI(10.0)<<"\t"<<Ion.getXI(11.0)<<"\t"<<Ion.getXI(12.0)<<"\t"<<Ion.getXI(13.0)<<"\t"<<Ion.getXI(14.0)<<"\t"<<Ion.getXI(15.0)<<"\t"<<Ion.getXI(16.0)<<"\t"<<Ion.getXI(17.0)<<"\t"<<Ion.getXI(18.0)<<"\t"<<Ion.getXI(19.0)<<"\t"<<Ion.getXI(20.0)<<endl;

		fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;
		z=6.0;
		while(z<20.2){
			fout<<"\t"<<Ion.getXI(z);
			z+=0.25;
		}
		fout<<endl;

}
fout.close();

*/
////////////////////////////////////////////////////////////////////////
// Grid likelihood values - step function
///////////////////////////////////////////////////////////////////////
/*
string files;
double zeta0, zeta1,z0,deltaz;
double like,z;
double p1min, p1max, p1step;
double p2min, p2max, p2step;

vector<double> param;
Ionization Ion(&c,&a);
int k;
int n1(151);
int n2(151);

z0=13.0;
deltaz=2.0;

p1min=30.0;
p1max=75.0;
p2min=0.0;
p2max=300.0;

p1step=(p1max-p1min)/(n1-1);
p2step=(p2max-p2min)/(n2-1);

files="likelihoodZeta_twostep.dat";
fout.open(files.c_str());

for(i=1;i<=n1;i++){
	zeta0=p1min+p1step*(double)(i-1);
	for(j=1;j<=n2;j++){
		zeta1=p2min+p2step*(double)(j-1);

		param.clear();
		param.push_back(2.0);
		param.push_back(zeta0);
		param.push_back(zeta1);
		param.push_back(z0);
		param.push_back(deltaz);
		Ion.setZetaParam(0.0,param,1);
		like=Ion.likelihoodZeta(param);
		cout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<like<<endl;
		fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;
		z=6.0;
		while(z<20.2){
			fout<<"\t"<<Ion.getXI(z);
			z+=0.25;
		}
		fout<<endl;


	}
}
fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Grid likelihood values - polynomial in Nion
///////////////////////////////////////////////////////////////////////
/*
string files;
double zeta0, zeta1,z0,deltaz;
double like, z;
vector<double> param;
Ionization Ion(&c,&a);
int k;
double clump;
double p1min,p1max,p1step;
double p2min, p2max, p2step;


int n1(151);
int n2(151);
//n1=20;
//n2=20;
//n1=50;
//n2=50;

//p1min=0.825;
//p1max=1.3;

//bolton
p1min=0.75;
p1max=1.55;

//fg
//p1min=0.83;
//p1max=1.15;

//p2min=7.0;
//p2max=20.0;

p2min=6.5;
p2max=21.0;

//wmap3
p2min=6.5;
p2max=30.0;

p1step=(p1max-p1min)/(double)(n1-1);
p2step=(p2max-p2min)/(double)(n2-1);

deltaz=2.0;

files="likelihoodZeta_ndot.dat";
fout.open(files.c_str());

for(i=1;i<=n1;i++){
	zeta0=p1min+(double)(i-1)*p1step;
	for(j=1;j<=n2;j++){
		zeta1=p2min+(double)(j-1)*p2step;

 		param.clear();
 		param.push_back(3);
 		param.push_back(zeta0);
 		param.push_back(zeta1);
 		param.push_back(0.0);
 		param.push_back(0.0);
		Ion.setZetaParam(0.0,param,1);
		like=Ion.likelihoodZeta(param);
 		cout<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like<<endl;
 		fout<<param[0]<<"\t"<<param[1]<<"\t"<<param[2]<<"\t"<<param[3]<<"\t"<<param[4]<<"\t"<<like;
 		z=6.0;
 		while(z<20.2){
 			fout<<"\t"<<Ion.getXI(z);
 			z+=0.25;
 		}
 		fout<<endl; 

	} 
} 
fout.close();
*/
//////////////////////////////////////////////////////////////////////
// Test 2D spline function
/////////////////////////////////////////////////////////////////////
/*
double x,y;
double *x1v, *x2v, **z;
int mm, nn;
double dx(0.5);
double dy(0.5);

Spline2D SP;

mm=10;
nn=10;

x1v=dvector(1,mm);
x2v=dvector(1,nn);
z=dmatrix(1,mm,1,nn);

//initialise
for(i=1;i<=mm;i++){
	x1v[i]=(double)(i-1)*dx;
	for(j=1;j<=nn;j++){
		x2v[j]=(double)(j-1)*dy;
		z[i][j]=pow(x1v[i],2.0)+pow(x2v[j],2.0);
	}
}

SP.setSplineSP(mm,nn,x1v,x2v,z);

cout<<SP.returnValue(3.3,2.4)<<endl;

*/
//////////////////////////////////////////////////////////////////////
// Test 2D clumping spline
/////////////////////////////////////////////////////////////////////
/*
double z, xi;
Ionization Ion(&c,&a);

z=ZION_START-0.1;

Ion.setClumping();

cout<<"here"<<endl;
cout<<Ion.getClumping(6.5,0.5)<<endl;

file="clumping.dat";
fout.open(file);

while(z>ZION_END){
	xi=0.9999;
	while(xi>0.01){
		cout<<z<<"\t"<<xi<<"\t"<<Ion.getClumping(z,xi)<<endl;
		fout<<z<<"\t"<<xi<<"\t"<<Ion.getClumping(z,xi)<<endl;
		xi/=1.1;
	}
	z-=0.3;
}

fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Convert Gamma constraints into Nion constraints
///////////////////////////////////////////////////////////////////////
    /*
double gamma0, gammam, gammap;
double sigm, sigp;
double z;
double nion0, nionm, nionp;
char tag[50];
 double nionNorm(1e51);

file="gamma_fg08.dat";
//file="gamma_bh07.dat";
fin.open(file);
fin.getline(tag,50);

file="nion_fg08.dat";
//file="nion_bh07.dat";
fout.open(file);

while(!fin.eof()){
	fin>>z>>gamma0>>sigm>>sigp;
	cout<<z<<"\t"<<gamma0<<"\t"<<sigm<<"\t"<<sigp<<endl;
	gammam=gamma0+sigm;
	gammap=gamma0+sigp;

	nion0=lyf.nionFromGamma(z,gamma0,3.0,3.0);
	nionm=lyf.nionFromGamma(z,gammam,3.0,3.0);
	nionp=lyf.nionFromGamma(z,gammap,3.0,3.0);

	//fout<<z<<"\t"<<nion0<<"\t"<<nionm-nion0<<"\t"<<nionp-nion0<<endl;
	fout<<z<<"\t"<<log10(nion0)<<"\t"<<log10(nionm)<<"\t"<<log10(nionp)<<endl;
}

fin.close();
fout.close();

 lyf.setLLSFlag(0);

file="gamma_fg08.dat";
//file="gamma_bh07.dat";
fin.open(file);
fin.getline(tag,50);

file="nion_fg08_nolls.dat";
//file="nion_bh07_nolls.dat";
fout.open(file);

while(!fin.eof()){
	fin>>z>>gamma0>>sigm>>sigp;
	cout<<z<<"\t"<<gamma0<<"\t"<<sigm<<"\t"<<sigp<<endl;
	gammam=gamma0+sigm;
	gammap=gamma0+sigp;

	nion0=lyf.nionFromGamma(z,gamma0,3.0,3.0);
	nionm=lyf.nionFromGamma(z,gammam,3.0,3.0);
	nionp=lyf.nionFromGamma(z,gammap,3.0,3.0);

	//fout<<z<<"\t"<<nion0<<"\t"<<nionm-nion0<<"\t"<<nionp-nion0<<endl;
	fout<<z<<"\t"<<log10(nion0)<<"\t"<<log10(nionm)<<"\t"<<log10(nionp)<<endl;
}

fin.close();
fout.close();
    */
////////////////////////////////////////////////////////////////////////
// Conversition table between ionization fraction and ionization fluctuation
///////////////////////////////////////////////////////////////////////
    /*
//this is the important bit for calculating the conversion table
    string files;
    int nx,nz;
    double z,dx,xi;
    double zmax, zmin, dz;
    double xmin, xmax;
    double **ftable;
    double *xiv, *zv;
	 double pTB, k(0.1);

	 nz=25;
	 zmax=18.0;
	 zmin=6.0;
	 dz=(zmax-zmin)/(double)(nz-1);
	 cout<<"delta_z= "<<dz<<endl;

	 nx=20;
	 xmax=1.0;
	 xmin=0.0;
	 dx=(xmax-xmin)/(double)(nx-1);

	 ftable=dmatrix(1,nx,1,nz);
	 xiv=dvector(1,nx);
	 zv=dvector(1,nz);

	 for(i=1;i<=nx;i++) xiv[i]=xmin+(i-1)*dx;
	 for(i=1;i<=nz;i++) zv[i]=zmin+(i-1)*dz;


	 files="fluc_key_short.dat";
	 fout.open(files.c_str());
	 


	   fout<<0.0;
	   for(i=0;i<nx;i++) fout<<"\t"<<i*dx;
	   fout<<endl;


	 for(j=1;j<=nz;j++){
	   z=zmin+(j-1)*dz;
	   xidd(0.0,z,&c,1);
	   fout<<z;
	   for(i=1;i<=nx;i++){
	     xi=(i-1)*dx;
	     pTB=rz.getPowerFZH(z,k,xi);  //fix xi=fcoll*zeta
	     fout<<"\t"<<pTB;
	     cout<<"\t"<<pTB<<endl;
	     ftable[i][j]=pTB;  
	   }
	   fout<<endl;
	   z-=1.0;
	 }
	 fout.close();

	 //use 2D spline to produce full table
	 Spline2D fspline;
	 int nbins;

	 fspline.setSplineSP(nx,nz,xiv,zv,ftable);

	 files="fluc_table_short.dat";
	 fout.open(files.c_str());

       nbins=100;
       dz=0.25;

       xmin=0.0;
       xmax=1.0;
       dx=(xmax-xmin)/(nbins-1);

       z=6.25;   
       while(z<18.0){

       for(i=0;i<nbins;i++){
	 xi=xmin+i*dx;
	 fout<<z<<"\t"<<xi<<"\t"<<fspline.returnValue(xi,z)<<endl;
       }

       z+=dz;
    }

	 fout.close();

	 free_dmatrix(ftable,1,nx,1,nz);
	 free_dvector(xiv,1,nx);
	 free_dvector(zv,1,nz);
    */
////////////////////////////////////////////////////////////////////////
// Optical depth for hitting a Lya opaque cloud
///////////////////////////////////////////////////////////////////////
    /*
     double ncloud,z;

     z=20.0;
     z=17.0;
     ncloud=nm(coolMass(&c,z),z,&c);

     cout<<z<<"\t"<<ncloud<<"\t"<<coolMass(&c,z)<<endl;
    */
////////////////////////////////////////////////////////////////////////
// Test random number generator
///////////////////////////////////////////////////////////////////////
    /*
    extern long iidum;

    iidum= -1;
    i=1;
    while(i>0){
      //cout<<dran2(&iidum)<<endl;
      i++;
      cout<<i<<endl;
    }
    */
////////////////////////////////////////////////////////////////////////
// Calculate the collapse fraction for different miniumum masses
///////////////////////////////////////////////////////////////////////
    /*
    double mmin, fcoll, z;
    double mcool, dfdz;
    double dz(0.01);
    
    z=30.0;
    file="fcoll_test.dat";
    fout.open(file);
    while(z>3.0){

      mcool=coolMass(&c,z);

      mmin=mcool;
      fcoll=c.fColl(z,mmin);
      dfdz=(c.fColl(z-dz,mmin)-c.fColl(z+dz,mmin))/(2.0*dz);
      cout<<z<<"\t"<<mmin<<"\t"<<fcoll<<"\t"<<dfdz<<endl;
      fout<<z<<"\t"<<mmin<<"\t"<<fcoll<<"\t"<<dfdz<<endl;
      z-=0.5;
    }
    fout.close();
    */
////////////////////////////////////////////////////////////////////////
// Test luminosity function calculations - galaxies
///////////////////////////////////////////////////////////////////////
/*
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);
    
    double alpha, phi0, L0;
    double M0, m0;
    double z;
    double Lmin;
    double nuObs;
    double e25;

    cout<<"Yoshida 2006"<<endl;
    z=4.7;
    alpha=-2.31;
    phi0=pow(10.0,5.76)/1.0e9;
    L0=pow(10.0,10.79)*LSOLAR;
    nuObs=SPEEDOFLIGHT_CGS/(1500e-8);

    Lmin=Gal.lumFromAbsMag(-20.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    Lmin=Gal.lumFromAbsMag(-18.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;



    cout<<"Bouwens 2006"<<endl;
    z=6.0;
    alpha=-1.73;
    phi0=pow(10.0,6.31)/1.0e9;
    L0=pow(10.0,10.50)*LSOLAR;
    nuObs=SPEEDOFLIGHT_CGS/(1350e-8);

    Lmin=Gal.lumFromAbsMag(-20.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    Lmin=Gal.lumFromAbsMag(-18.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    cout<<"Bunker 2004"<<endl;
    z=6.0;
    alpha=-1.6;
    phi0=pow(10.0,5.36)/1.0e9;
    L0=pow(10.0,10.75)*LSOLAR;
    nuObs=SPEEDOFLIGHT_CGS/(1350e-8);

    Lmin=Gal.lumFromAbsMag(-20.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    Lmin=Gal.lumFromAbsMag(-18.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;


    cout<<"Other data sets"<<endl;
    cout<<"Steidel 1999"<<endl;
    z=3.0;
    alpha=-1.6;
    phi0=1.6e-2;
    m0=24.48;
    nuObs=SPEEDOFLIGHT_CGS/(1350e-8);
    L0=Gal.lumFromMag(m0,z)*nuObs;
    Lmin=0.1*L0;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    cout<<"Bouwens 2009"<<endl;
    z=3.8;
    alpha=-1.73;
    phi0=1.3e-3;
    M0= -20.98;
    nuObs=SPEEDOFLIGHT_CGS/(1600e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;
    Lmin=Gal.lumFromAbsMag(-16.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    z=5.0;
    alpha=-1.66;
    phi0=1.0e-3;
    M0= -20.64;
    nuObs=SPEEDOFLIGHT_CGS/(1600e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;
    Lmin=Gal.lumFromAbsMag(-16.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    z=5.9;
    alpha=-1.74;
    phi0=1.4e-3;
    M0= -20.24;
    nuObs=SPEEDOFLIGHT_CGS/(1350e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;
    Lmin=Gal.lumFromAbsMag(-16.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    z=7.4;
    alpha=-1.74;
    phi0=1.4e-3;
    M0= -19.3;
    nuObs=SPEEDOFLIGHT_CGS/(1900.0e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;
    Lmin=Gal.lumFromAbsMag(-16.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;

    //Bouwyns 0912.4263
    z=10.2;
    alpha=-1.74;
    phi0=1.2e-3;
    M0= -18.7;
    nuObs=SPEEDOFLIGHT_CGS/(1600.0e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;
    Lmin=Gal.lumFromAbsMag(-16.0,z)*nuObs;
    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;
*/
////////////////////////////////////////////////////////////////////////
// Test luminosity function calculations - lower limit of integration
///////////////////////////////////////////////////////////////////////
/*
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);
    stringstream ss;
    string files;
    double zin, lamObs;
    double x;
    
    double alpha, phi0, L0;
    double M0, m0;
    double z;
    double Lmin;
    double nuObs;
    double e25;
    double sig_alpha;
    double sig_M0;
    double sig_phi0;

    z=zin_arg;

     //redshift dependent astrophysics parameters
      ss.clear();
       ss<<"./lf_z"<<z<<".dat";
       ss>>files;
       cout<<files<<endl;
       fin.open(files.c_str());
       fin>>zin>>lamObs;
       fin>>alpha>>sig_alpha;
       fin>>phi0>>sig_phi0;
       fin>>M0>>sig_M0;
       fin.close();
       z=zin;

   nuObs=SPEEDOFLIGHT_CGS/(1900.0e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;

    x=10.0;

    ss.clear();
    ss<<"./e25_lmin_z"<<z<<".dat";
    ss>>files;
    fout.open(files.c_str());

       while(x>1.0e-4){
          Lmin=x*L0;

    e25=Gal.emissivityGal(Lmin,nuObs,z,phi0,L0,alpha);
    fout<<z<<"\t"<<x<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<Lmin<<"\t"<<e25<<endl;
    x/=2.0;
       }
       fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Test luminosity function calculations -QSO
///////////////////////////////////////////////////////////////////////
/*
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);
    
    double beta1, beta2, phi0, L0;
    double M0, m0;
    double z;
    double Lmin;
    double nuObs;
    double e25;

    cout<<"QS"<<endl;
    z=5.0;
    beta1= -1.24;
    beta2= -2.70;
    //beta2= -2.20;
    phi0=pow(10.0,2.51)/1.0e9;
    L0=pow(10.0,11.89)*LSOLAR;
    nuObs=SPEEDOFLIGHT_CGS/(1450e-8);

    Lmin=Gal.lumFromAbsMag(-22.0,z)*nuObs;
    e25=Gal.emissivityQSO(Lmin,nuObs,z,phi0,L0,beta1,beta2);
    cout<<z<<"\t"<<nuObs<<"\t"<<phi0<<"\t"<<log10(Lmin/LSOLAR)<<"\t"<<Lmin/L0<<"\t"<<e25<<endl;

    z=6.0;
    beta1= -1.24;
    beta2= -2.70;
    //beta2= -2.20;
    phi0=pow(10.0,2.51)/1.0e9;
    L0=pow(10.0,11.53)*LSOLAR;
    nuObs=SPEEDOFLIGHT_CGS/(1450e-8);

    cout<<Gal.lumFromAbsMag(-21.4,z)*nuObs/L0<<endl;

    Lmin=Gal.lumFromAbsMag(-22.0,z)*nuObs;
    e25=Gal.emissivityQSO(Lmin,nuObs,z,phi0,L0,beta1,beta2);
    cout<<z<<"\t"<<nuObs<<"\t"<<L0<<"\t"<<phi0<<"\t"<<log10(Lmin/LSOLAR)<<"\t"<<Lmin/L0<<"\t"<<e25<<endl;

    double gamma1, gamma2;
    cout<<"Hopkins"<<endl;
    z=5.0;
    gamma1=0.497;
    gamma2=1.57;
    beta1= -(gamma1+1.0);
    beta2= -(gamma2+1.0);
    phi0=pow(10.0,-5.38)/log(10.0);
    L0=pow(10.0,12.49)*LSOLAR;
    nuObs=SPEEDOFLIGHT_CGS/(1450e-8);

    cout<<Gal.lumFromAbsMag(-21.4,z)*nuObs/L0<<endl;

    Lmin=Gal.lumFromAbsMag(-22.0,z)*nuObs;
    e25=Gal.emissivityQSO(Lmin,nuObs,z,phi0,L0,beta1,beta2);
    cout<<z<<"\t"<<nuObs<<"\t"<<L0<<"\t"<<phi0<<"\t"<<log10(Lmin/LSOLAR)<<"\t"<<Lmin/L0<<"\t"<<beta1<<"\t"<<beta2<<"\t"<<e25<<endl;

*/

////////////////////////////////////////////////////////////////////////
// Test Stark+ LF calculation
///////////////////////////////////////////////////////////////////////
/*
     double z, L, nL, mag, sfr, nLdO, dV;
     double magmin, magmax, mstep;
     double magapp;
     double mass, massC, massI;
     int nstep(20);
     int nz(10);
     string files;
     stringstream ss;
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);

   
    for(j=0;j<=nz;j++){
    z=4.0+j*1.0;

    magmin= -14.0;
    magmax= -28.0;

    ss.clear();
    ss<<"./stark_lf_z"<<z<<".dat";
    ss>>files;
    fout.open(files.c_str());

    mstep=(magmax-magmin)/(double)(nstep-1);
    for(i=0;i<=nstep;i++){
       mag=magmin+i*mstep;
       L=Gal.lumFromAbsMag(mag,z);
       magapp=Gal.magFromL(L,z);
       sfr=Gal.sfrFromLum(L);
       nL=Gal.starkLF(L,z);
       dV=c.volumeComoving(z);
       nLdO=nL*dV*pow(PI/180.0,2.0);
       mass=Gal.massFromLum(L,z);
       massC=coolMass(&c,z);
       massI=minIonMass(&c,z);
       fout<<mag<<"\t"<<L<<"\t"<<sfr<<"\t"<<nL<<"\t"<<nLdO<<"\t"<<magapp<<"\t"<<mass<<"\t"<<mass/massC<<"\t"<<mass/massI<<endl;
    }
    fout.close();
    
    }
*/
////////////////////////////////////////////////////////////////////////
// Calculate recombination rate for maintaining reionization
///////////////////////////////////////////////////////////////////////
/*
     double z, recomb, C;
     double ndot, alphaS(3.0);
     double e25crit;
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);

     z=6.0;
     fesc=1.0;
     C=1;

     z=3.8;
     ndot=Gal.ndotGal(1.0,alphaS,fesc);
     recomb=Gal.recombLevel(z,1.0,C);
     e25crit=recomb/ndot;
     cout<<z<<"\t"<<recomb<<"\t"<<ndot<<"\t"<<e25crit<<endl;

     z=5.0;
     ndot=Gal.ndotGal(1.0,alphaS,fesc);
     recomb=Gal.recombLevel(z,1.0,C);
     e25crit=recomb/ndot;
     cout<<z<<"\t"<<recomb<<"\t"<<ndot<<"\t"<<e25crit<<endl;
     
     z=5.9;
     ndot=Gal.ndotGal(1.0,alphaS,fesc);
     recomb=Gal.recombLevel(z,1.0,C);
     e25crit=recomb/ndot;
     cout<<z<<"\t"<<recomb<<"\t"<<ndot<<"\t"<<e25crit<<endl;
  
     z=7.4;
     ndot=Gal.ndotGal(1.0,alphaS,fesc);
     recomb=Gal.recombLevel(z,1.0,C);
     e25crit=recomb/ndot;
     cout<<z<<"\t"<<recomb<<"\t"<<ndot<<"\t"<<e25crit<<endl;

     z=9.0;
     ndot=Gal.ndotGal(1.0,alphaS,fesc);
     recomb=Gal.recombLevel(z,1.0,C);
     e25crit=recomb/ndot;
     cout<<z<<"\t"<<recomb<<"\t"<<ndot<<"\t"<<e25crit<<endl;
*/
////////////////////////////////////////////////////////////////////////
// Compare galaxy emissivity and recombination barrier assuming Stark+ model
///////////////////////////////////////////////////////////////////////
/*
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);
    stringstream ss;
    string files;
    double zin, lamObs;
    double x;
    double mmin;
    
    double alpha, phi0, L0;
    double M0, m0;
    double z;
    double Lmin, LminC, LminI, LminM;
    double ndotC, ndotI, ndotM;
    double nuObs;
    double e25, e25I, e25C, e25M;
    double sig_alpha;
    double sig_M0;
    double sig_phi0;
    double ndot, recomb, e25crit;
    double C(1.0), alphaS(3.0);
    fesc=1.0;
    vector<double> lfparam;

    files="./recomb_lim.dat";
    fout.open(files.c_str());

    z=5.0;
    while(z<16.0){


    lfparam=Gal.getLFParam(z);
    M0=lfparam[0];
    phi0=lfparam[1];
    alpha=lfparam[2];
    nuObs=SPEEDOFLIGHT_CGS/(1900.0e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;
    cout<<L0<<endl;



    mmin=coolMass(&c,z);
    LminC=Gal.lumFromSFR(Gal.sfrHalo(mmin,z))*nuObs;
    e25C=Gal.emissivityGal(LminC,nuObs,z,phi0,L0,alpha);

    mmin=minIonMass(&c,z);
    LminI=Gal.lumFromSFR(Gal.sfrHalo(mmin,z))*nuObs;
    e25I=Gal.emissivityGal(LminI,nuObs,z,phi0,L0,alpha);

    mmin=coolMassH2(&c,z);
    LminM=Gal.lumFromSFR(Gal.sfrHalo(mmin,z))*nuObs;
    e25M=Gal.emissivityGal(LminM,nuObs,z,phi0,L0,alpha);


     ndotC=Gal.ndotGal(e25C,alphaS,fesc);
     ndotI=Gal.ndotGal(e25I,alphaS,fesc);
     ndotM=Gal.ndotGal(e25M,alphaS,fesc);
     ndot=Gal.ndotGal(1.0,alphaS,fesc);
     recomb=Gal.recombLevel(z,1.0,C);
     e25crit=recomb/ndot;

     cout<<z<<"\t"<<LminC<<"\t"<<LminI<<"\t"<<LminM<<"\t"<<L0<<"\t"<<e25C<<"\t"<<e25I<<"\t"<<e25M<<"\t"<<e25crit<<"\t"<<ndotC<<"\t"<<ndotI<<"\t"<<ndotM<<"\t"<<recomb<<endl;
     fout<<z<<"\t"<<LminC<<"\t"<<LminI<<"\t"<<LminM<<"\t"<<L0<<"\t"<<e25C<<"\t"<<e25I<<"\t"<<e25M<<"\t"<<e25crit<<"\t"<<ndotC<<"\t"<<ndotI<<"\t"<<ndotM<<"\t"<<recomb<<endl;
     z+=1.0;
    }
    fout.close();
*/

////////////////////////////////////////////////////////////////////////
// Calculate Stellar density needed to ionize universe at a given redshift
///////////////////////////////////////////////////////////////////////
/*
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);
    stringstream ss;
    string files;
    double sfd, z, sfr, sfrZ;
    double dzdt;
     double sfdR;
     double dz(0.5);

     file="./sfd_reion.dat";
     fout.open(file);

     sfdR=Gal.reionSFD(1.0);

     z=4.0;
     while(z<31.0){
        sfr=Gal.sfrDensity(z);
        sfd=Gal.totalSFD(z);
        dzdt=(1.0+z)*c.hubbleZ(z)*c.getH()*UNH*YEAR;
        sfrZ=sfr/dzdt;
        cout<<z<<"\t"<<sfr<<"\t"<<sfd<<"\t"<<sfdR<<"\t"<<sfrZ<<endl;
        fout<<z<<"\t"<<sfr<<"\t"<<sfd<<"\t"<<sfdR<<"\t"<<sfrZ<<endl;
        z+=dz;
     }
     fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Check halo abundance available for forming galaxies
///////////////////////////////////////////////////////////////////////
/*
     double nHI, nH2;
     double mHI, mH2;
     double mIon, nIon;
     double ngal;
     double z;
    double mmin;
    
    double alpha, phi0, L0;
    double M0, m0;
    double nuObs, LminC;
    vector<double> lfparam;
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);

     files="./halo_abundance.dat";
     fout.open(files.c_str());

     z=40.0;
     while(z>3.0){
        mHI=coolMass(&c,z);
        mH2=coolMassH2(&c,z);
        mIon=minIonMass(&c,z);
        nHI=c.nCollObject(z,mHI);
        nH2=c.nCollObject(z,mH2);
        nIon=c.nCollObject(z,mIon);

        //now galaxies from LF
    lfparam=Gal.getLFParam(z);
    M0=lfparam[0];
    phi0=lfparam[1];
    alpha=lfparam[2];
    nuObs=SPEEDOFLIGHT_CGS/(1500.0e-8);
    L0=Gal.lumFromAbsMag(M0,z)*nuObs;
    mmin=coolMass(&c,z);
    LminC=Gal.lumFromSFR(Gal.sfrHalo(mmin,z))*nuObs;
    //cout<<LminC<<"\t"<<L0<<"\t"<<LminC/L0<<endl;
    //cout<<phi0<<"\t"<<alpha<<endl;
    ngal=Gal.numberDensityLF(LminC,alpha,phi0,L0);

    cout<<z<<"\t"<<mHI<<"\t"<<mH2<<"\t"<<mIon<<"\t"<<nHI<<"\t"<<nH2<<"\t"<<nIon<<"\t"<<ngal<<"\t"<<phi0<<endl;
    fout<<z<<"\t"<<mHI<<"\t"<<mH2<<"\t"<<mIon<<"\t"<<nHI<<"\t"<<nH2<<"\t"<<nIon<<"\t"<<ngal<<"\t"<<phi0<<endl;
    z-=0.25;
     }


     fout.close();
*/
////////////////////////////////////////////////////////////////////////
// Radio loud number of quasars
///////////////////////////////////////////////////////////////////////
/*
     double nqso;
     double z;
    
    double alpha, phi0, L0;
    double M0, m0;
    double nuObs, LminC;
    vector<double> lfparam;
    Galaxies Gal;
    Gal.initGalaxies(&c,&a);

     files="./halo_abundance.dat";
     fout.open(files.c_str());

     nqso=Gal.luminosityFunctionQSO_Hop(L,phi0,L0,gamma1,gamma2);


     fout.close();

*/
////////////////////////////////////////////////////////////////////////
// testing pdf code
///////////////////////////////////////////////////////////////////////
/*
     double mpix(1.0e13);
     double z, zeta, psi;
     double xi, delta, pdf;
     z=20.0;
     zeta=40.0;
     psi=1.0;
     delta=0.0;

     pdf=tocm.getPsiPDF(psi,z,zeta,mpix);
     xi=tocm.XIfromOverdensity(delta,z,zeta,mpix);

     delta= -1.0;
     file="./psi_temp.dat";
     fout.open(file);
     while(delta<1.4){
        psi=tocm.psiFromDelta(delta,z,zeta,mpix);
        xi=tocm.XIfromOverdensity(delta,z,zeta,mpix);
        cout<<delta<<"\t"<<psi<<"\t"<<xi<<endl;
        fout<<delta<<"\t"<<psi<<"\t"<<xi<<"\t"<<delta/c.growthFac(z)<<endl;
        delta+=0.025;
     }
     fout.close();

     psi=1.;
     tocm.deltaFromPsi(psi,z,zeta,mpix);

     psi=0.01;
     file="./pdf_temp.dat";
     fout.open(file);
     while(psi<1.4){
        pdf=tocm.getBrightnessPDF(psi,z,zeta,mpix);
        cout<<psi<<"\t"<<pdf<<endl;
        fout<<psi<<"\t"<<pdf<<endl;
        psi+=0.05;
     }
     fout.close();

     cout<<xi<<endl;

     double prob;
     prob=tocm.evaluateCondProb(delta,z,zeta,mpix);
     cout<<prob<<endl;
*/
//////////////////////////////////////////////////////////////////
// Calculate the number of galaxies within an observed beam
//////////////////////////////////////////////////////////////////
/*
     double z, mass, N,N2;

     z=8.0;
     mass=coolMass(&c,z);
     //N=meanHaloNumber(z,&c,mass);
     N=c.nCollObject(z);
     N2=c.nCollObject(z,1.0e10);
     cout<<z<<"\t"<<mass<<"\t"<<N<<"\t"<<N2<<endl;
*/
//////////////////////////////////////////////////////////////////
// Test bolometric to band conversions from Hopkins (2007)
//////////////////////////////////////////////////////////////////
/*
  Galaxies Gal;
  stringstream ss;
  FisherLF lf;
  CosmParam fiducial;
  string namec="_fiducial";

  file="./project.ini";
  fiducial=lf.setCosmParamFromFile(file,namec);
  Cosmology cosm(fiducial.omm,fiducial.oml,fiducial.omb,fiducial.h,fiducial.sigma8,fiducial.nscal,fiducial.omnu);
  Gal.initGalaxies(&cosm,NULL);


  double Lmax(14.0), Lmin(8.0);
  double Lstep;
  int nL(30);
  double L, Lband, Ltest;
  double nu(1.0);

  Lstep=(Lmax-Lmin)/(double)(nL);

  L=Lmin;
  files="bolo_test.dat";
  fout.open(files.c_str());

  for(i=0;i<=nL;i++){
     L=Lmin+i*Lstep;
     L=pow(10.0,L)*LSOLAR;

     Lband=Gal.boloToBandL(nu,L);
     Ltest=Gal.bandToBoloL(nu,Lband);
     fout<<L/LSOLAR<<"\t"<<Lband/LSOLAR<<"\t"<<Ltest/LSOLAR<<endl;

  }
  fout.close();
*/
//////////////////////////////////////////////////////////////////
// Test QSO LF parameter fits from Hopkins (2007)
//////////////////////////////////////////////////////////////////
/*
  Galaxies Gal;
  stringstream ss;
  FisherLF lf;
  CosmParam fiducial;
  string namec="_fiducial";

  file="./project.ini";
  fiducial=lf.setCosmParamFromFile(file,namec);
  Cosmology cosm(fiducial.omm,fiducial.oml,fiducial.omb,fiducial.h,fiducial.sigma8,fiducial.nscal,fiducial.omnu);
  Gal.initGalaxies(&cosm,NULL);

  vector<double> LFparam;
  double lphi0, lL0, gamma1, gamma2; //QSO LF parameters
  double ndens, Lmin, Mmin;
  double Lsol(3.9e33);
  double z, zmin,zmax;
  double dz(0.5);
  int nz(40);
  double nu(1.0),LBmin;

  zmin=0.1;
  zmax=6.0;
  dz=(zmax-zmin)/(double)(nz);

  //Convert from a B-band magnitude to the QSO luminosity
  //not a simple conversion from Abs Mag to L as code currently does!!!!
  Mmin= -27.0;

  files="./qso_lf_z_ple.dat";
  fout.open(files.c_str());
  for(i=0;i<=nz;i++){
     z=zmin+i*dz;
     LFparam=Gal.getLFParamQSO_PLE(z);   
     lphi0=LFparam[0];
     lL0=LFparam[1];
     gamma1=LFparam[2];
     gamma2=LFparam[3];
     //LBmin=Gal.lumFromAbsMag(Mmin,z); //get B band L
     //Lmin=Gal.bandToBoloL(nu,LBmin);
     Lmin=1.1e13;

     //cout<<LBmin/LSOLAR<<"\t"<<Lmin/LSOLAR<<endl;
     ndens=Gal.numberDensityQSO_HOP(z,Lmin,1);
     ndens=log10(ndens);
     cout<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<ndens<<endl;
     fout<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<ndens<<endl;
  }
  fout.close();

  files="./qso_lf_z_full.dat";
  fout.open(files.c_str());
  for(i=0;i<=nz;i++){
     z=zmin+i*dz;
     LFparam=Gal.getLFParamQSO_FULL(z);   
     lphi0=LFparam[0];
     lL0=LFparam[1];
     gamma1=LFparam[2];
     gamma2=LFparam[3];
     //LBmin=Gal.lumFromAbsMag(Mmin,z); //get B band L
     //Lmin=Gal.bandToBoloL(nu,LBmin);
     Lmin=1.0e13;
     //Lmin=1.0e14; //quite close to Hopkins plot

     //cout<<Lmin<<"\t"<<pow(10.0,lL0)<<endl;
     ndens=Gal.numberDensityQSO_HOP(z,Lmin,0);
     ndens=log10(ndens);
     cout<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<ndens<<endl;
     fout<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<ndens<<endl;
  }
  fout.close();

*/
//////////////////////////////////////////////////////////////////
// Plot QSO LF  from Hopkins (2007)
// i.e. \phi(L) against L
//////////////////////////////////////////////////////////////////
/*
// THIS IS THE ONE TO UNCOMMENT!!!!!

  Galaxies Gal;
  stringstream ss;
  FisherLF lf;
  CosmParam fiducial;
  string namec="_fiducial";

     Astronomy Ast;
     MAGBANDS band;
     double m, flux;
     double dLmag, cH, Lnu, zmag;
     double nu, Lbol, Lband, lambda, dnu, magB;
     band=B;
     zmag=0.0;
     dLmag=10.0/1.0e6; //10 pc, so get absolute mag!
     m= -27.0;
     nu=SPEEDOFLIGHT_CGS/4400.0e-10;





  file="./project.ini";
  fiducial=lf.setCosmParamFromFile(file,namec);
  Cosmology cosm(fiducial.omm,fiducial.oml,fiducial.omb,fiducial.h,fiducial.sigma8,fiducial.nscal,fiducial.omnu);
  Gal.initGalaxies(&cosm,NULL);

  vector<double> LFparam;
  double lphi0, lL0, gamma1, gamma2; //QSO LF parameters
  double ndens, Lmin, Mmin, Lmax;
  double Lsol(3.9e33);
  double z, zmin,zmax;
  double dz(0.5);
  int nz(5);
  //double nu(1.0),LBmin;
  int nL(20);
  double lLmin, lLmax,lstep;
  double L,dndL,dndLB,phi0,L0;
  zmin=1.0;
  zmax=6.0;
  dz=1.0;

  //Convert from a B-band magnitude to the QSO luminosity
  //not a simple conversion from Abs Mag to L as code currently does!!!!
  Mmin= -27.0;

  for(i=0;i<=nz;i++){
     z=zmin+i*dz;
     files="./qso_lf_L_ple_z";
     ss.clear();
     ss<<files<<z<<".dat";
     ss>>files;
     fout.open(files.c_str());

     LFparam=Gal.getLFParamQSO_PLE(z);   
     lphi0=LFparam[0];
     lL0=LFparam[1];
     gamma1=LFparam[2];
     gamma2=LFparam[3];
     Lmin=1.0e8;
     Lmax=1.0e17;

     lLmin=log10(Lmin);
     lLmax=log10(Lmax);
     lstep=(lLmax-lLmin)/(nL);

     phi0=pow(10.0,lphi0);
     L0=pow(10.0,lL0);

     cout<<"# "<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<phi0<<"\t"<<L0<<endl;
     fout<<"# "<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<phi0<<"\t"<<L0<<endl;
     for(int j=0;j<=nL;j++){
        L=pow(10.0,lLmin+lstep*j);
        dndL=Gal.luminosityFunctionQSO_dNdLogL(L,phi0,L0,gamma1,gamma2);

        Lband=Gal.boloToBandL(nu,L*LSOLAR)*Gal.attenuationFactor(nu,L*LSOLAR); 
        lambda=effwave[band]*1.0e-8;
        dnu=deltawave[band]*1.0e-8*SPEEDOFLIGHT_CGS/lambda/lambda;
        //dnu=SPEEDOFLIGHT_CGS/(deltawave[band]*1.0e-8);
        Lnu=Lband/dnu;
        magB=Ast.magFromL(Lnu);
        dndLB=Gal.boloToBandPhi(nu,L*LSOLAR,dndL);

        cout<<L<<"\t"<<dndL<<"\t"<<log10(dndL)<<"\t"<<magB<<"\t"<<dndLB<<endl;
        fout<<L<<"\t"<<dndL<<"\t"<<log10(dndL)<<"\t"<<magB<<"\t"<<dndLB<<endl;
     }
     fout.close();
  }


  for(i=0;i<=nz;i++){
     z=zmin+i*dz;
     files="./qso_lf_L_full_z";
     ss.clear();
     ss<<files<<z<<".dat";
     ss>>files;
     fout.open(files.c_str());

     LFparam=Gal.getLFParamQSO_FULL(z);   
     lphi0=LFparam[0];
     lL0=LFparam[1];
     gamma1=LFparam[2];
     gamma2=LFparam[3];
     Lmin=1.0e8;
     Lmax=1.0e17;

     lLmin=log10(Lmin);
     lLmax=log10(Lmax);
     lstep=(lLmax-lLmin)/(nL);

     phi0=pow(10.0,lphi0);
     L0=pow(10.0,lL0);

     cout<<"# "<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<phi0<<"\t"<<L0<<endl;
     fout<<"# "<<z<<"\t"<<lL0<<"\t"<<lphi0<<"\t"<<gamma1<<"\t"<<gamma2<<"\t"<<phi0<<"\t"<<L0<<endl;
     for(int j=0;j<=nL;j++){
        L=pow(10.0,lLmin+lstep*j);
        dndL=Gal.luminosityFunctionQSO_dNdLogL(L,phi0,L0,gamma1,gamma2);

        Lband=Gal.boloToBandL(nu,L*LSOLAR)*Gal.attenuationFactor(nu,L*LSOLAR); 
        lambda=effwave[band]*1.0e-8;
        dnu=deltawave[band]*1.0e-8*SPEEDOFLIGHT_CGS/lambda/lambda;
        //dnu=SPEEDOFLIGHT_CGS/(deltawave[band]*1.0e-8);
        Lnu=Lband/dnu;
        magB=Ast.magFromL(Lnu);
        cout<<L*LSOLAR<<"\t"<<Lband<<"\t"<<Lnu<<"\t"<<dnu<<"\t"<<lambda<<endl;

        dndLB=Gal.boloToBandPhi(nu,L*LSOLAR,dndL);

        cout<<L<<"\t"<<dndL<<"\t"<<log10(dndL)<<"\t"<<magB<<"\t"<<dndLB<<endl;
        fout<<L<<"\t"<<dndL<<"\t"<<log10(dndL)<<"\t"<<magB<<"\t"<<dndLB<<endl;
     }
     fout.close();
  }


*/
///////////////////////////////////////////////////////////////////////
// Test astronomy conversions for different optical bands
////////////////////////////////////////////////////////////////////////
/*
     Astronomy Ast;
     MAGBANDS band;
     double m, flux;

     //
  Galaxies Gal;
  stringstream ss;
  FisherLF lf;
  CosmParam fiducial;
  string namec="_fiducial";

  file="./project.ini";
  fiducial=lf.setCosmParamFromFile(file,namec);
  Cosmology cosm(fiducial.omm,fiducial.oml,fiducial.omb,fiducial.h,fiducial.sigma8,fiducial.nscal,fiducial.omnu);
  Gal.initGalaxies(&cosm,NULL);
     //

     m=0.0;

     band=V;
     flux=Ast.fluxFromMagnitude(m,band);
     cout<<band<<"\t"<<m<<"\t"<<flux<<endl;

     for(i=0;i<=8;i++){
        band=MAGBANDS(i);
     flux=Ast.fluxFromMagnitude(m,band);
     cout<<band<<"\t"<<m<<"\t"<<flux<<endl;
     }

     //get the luminosity correponding to given magnitude
     double dL, cH, Lband, z;
     double nu, Lbol;
     band=B;
     z=0.0;
     dL=10.0/1.0e6; //10 pc, so get absolute mag!

     m= -27.0;
     Lband=Ast.lumBandFromMagnitude(m,band,dL,z);
     nu=SPEEDOFLIGHT_CGS/4400.0e-10;
     Lbol=Gal.bandToBoloL(nu,Lband);
     cout<<"Lband="<<Lband<<endl;
     cout<<"Lbol="<<Lbol<<"\t"<<Lbol/LSOLAR<<endl;    

     files="./Bmagnitude.dat";
     fout.open(files.c_str());
     Lbol=0.0;
     m= -12.0;
     while(m> -40.0){
        Lband=Ast.lumBandFromMagnitude(m,band,dL,z);
        cout<<m<<"\t"<<Lband<<endl;
        nu=SPEEDOFLIGHT_CGS/4400.0e-10;

        //correct for attenuation

        Lbol=Gal.bandToBoloL(nu,Lband);
        cout<<"Lband="<<Lband<<endl;
        cout<<"Lbol="<<Lbol<<"\t"<<Lbol/LSOLAR<<endl; 
        fout<<m<<"\t"<<Lband/LSOLAR<<"\t"<<Lbol/LSOLAR<<"\t"<<Lbol/Lband<<endl;
        m-=0.5;
     }
     fout.close();
*/
///////////////////////////////////////////////////////////////////////
  return 0;
}




