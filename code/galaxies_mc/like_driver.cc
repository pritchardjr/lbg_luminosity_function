
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
     int nthread(0);

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
      lls_flag=atoi(&argv[1][2]);
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
int ndot_flag(0), cmb_flag(0);
int param_flag;

param_flag=3;  //Nion polynomial
ndot_flag=0;
cmb_flag=0;

int n1(151);
int n2(151);
//n1=20;
//n2=20;
//n1=50;
//n2=50;

//ndot

//bolton
p1min=0.75;
p1max=1.55;

if(ndot_flag==1){
//fg
//p1min=0.83;
//p1max=1.15;
}

//cmb

p2min=6.5;
p2max=21.0;

//wmap3
if(cmb_flag==1){
p2min=6.5;
p2max=30.0;
}

p1step=(p1max-p1min)/(double)(n1-1);
p2step=(p2max-p2min)/(double)(n2-1);

deltaz=2.0;

Ion.setFlags(ndot_flag,cmb_flag);

files="likelihoodZeta_ndot.dat";
fout.open(files.c_str());

for(i=1;i<=n1;i++){
	zeta0=p1min+(double)(i-1)*p1step;
	for(j=1;j<=n2;j++){
		zeta1=p2min+(double)(j-1)*p2step;

 		param.clear();
 		param.push_back(param_flag);
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
////////////////////////////////////////////////////////////////////////
// Grid likelihood values - general case
///////////////////////////////////////////////////////////////////////
     string data_dir;
     string file_like;
     string files;
     double zeta0, zeta1,z0,deltaz;
     double zeta2, zeta3;
     double like, z;
     vector<double> param;
     Ionization Ion(&c,&a);
     //int k;
     //double clump;
     double p1min(0.0),p1max(0.0),p1step;
     double p2min(0.0), p2max(0.0), p2step;     

     int n1(151);
     int n2(151);
     n1=20;
     n2=20;
//n1=50;
//n2=50;
// n1=5;
// n2=5;

     //param_flag 1: ??   2: ??  3: Nion
   //cmb_flag  0 WMAP5  1 WMAP3  2 PLANCK
   //ndot_flag 0 BOLTON 1 FG

     //param_flag=3;  //Nion polynomial
     //ndot_flag=0;
     //cmb_flag=0;

     data_dir="./";
    file_like="likelihoodZeta.dat";

    if(param_flag==2) data_dir=data_dir+"twostep/";
    if(param_flag==3) data_dir=data_dir+"ndot/";

    if(cmb_flag==1) data_dir=data_dir+"wmap3/";
    if(cmb_flag==2) data_dir=data_dir+"planck/";
    if(ndot_flag==0 && cmb_flag==0) data_dir=data_dir+"bolton/"; 
    if(ndot_flag==1) data_dir=data_dir+"fg_data/"; 

    if(lls_flag==0) data_dir=data_dir+"no_lls/";

    files="mkdir "+data_dir;
    system(files.c_str());

    if(param_flag==2){
       z0=13.0;
       deltaz=2.0;

       p1min=30.0;
       p1max=75.0;
       p2min=0.0;
       p2max=300.0;

       zeta2=z0;
       zeta3=deltaz;
    }

     if(param_flag==3){  //Nion polynomial
        //ndot

        //bolton
        p1min=0.75;
        p1max=1.55;

        if(ndot_flag==1){
        //fg
           p1min=0.83;
           p1max=1.15;
        }

        //cmb

        p2min=6.5;
        p2max=21.0;

        //wmap3
        if(cmb_flag==1){
           p2min=6.5;
           p2max=30.0;
        }

	if(lls_flag==0){
	  p1min/=3.0;
	  p2max/=3.0;
	}

        zeta2=0.0;
        zeta3=0.0;
     }

p1step=(p1max-p1min)/(double)(n1-1);
p2step=(p2max-p2min)/(double)(n2-1);

deltaz=2.0;

 Ion.setFlags(ndot_flag,cmb_flag,lls_flag,nthread);
 lyf.setLLSFlag(lls_flag);



 files=data_dir+file_like;
 cout<<files<<endl;
fout.open(files.c_str());

for(i=1;i<=n1;i++){
	zeta0=p1min+(double)(i-1)*p1step;
	for(j=1;j<=n2;j++){
		zeta1=p2min+(double)(j-1)*p2step;

 		param.clear();
 		param.push_back(param_flag);
 		param.push_back(zeta0);
 		param.push_back(zeta1);
 		param.push_back(zeta2);
 		param.push_back(zeta3);
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
 


////////////////////////////////////////////////////////////////////////
// Process likelihood data to get neutral fraction distribution
///////////////////////////////////////////////////////////////////////

   //string data_dir;
   // string file_like;
     //string files;
     string str;
     int npoint(0);
     int ndata(0);
     int nentry;
     int nbins,nz,ibin,np1,np2;
     int iz;
     double pmin, pmax;
     double p1, p2, dp1, dp2, trash;
     //double like;
     vector<double> p1v, p2v,likev, xv;
     vector<double> pdf;
     //double p1min(1e30), p1max(-1e30), p2min(1e30), p2max(-1e30);
     //char line[5];
     double likeTot(0.0);
     double xmin, xmax, dx;
     double xi,xiold,xinew;
     //double z;
     double wnorm;
     double zmin(6.0);
     double dz(0.25);
     double likeMax(0.0);
     double p3, p4, p0;
     double p1Best, p2Best,p3Best,p4Best;

    p1min=1e30;
    p1max= -1e30;
    p2min=1e30;
    p2max= -1e30;


//marginalised parameter pdfs

    files=data_dir+file_like;
    cout<<files<<endl;

    fin.open(files.c_str()); 
    while(!fin.eof()){
      fin>>pmin;
       npoint++;
    }
    fin.close();
    npoint--;

    cout<<npoint<<endl;

    fin.open(files.c_str()); 
    while(!fin.eof()){
       getline(fin,str);
       ndata++;
    }
    fin.close();
    ndata--;

    cout<<npoint<<"\t"<<ndata<<endl;
    //find max and min values of parameters and likelihood total
    fin.open(files.c_str()); 
    while(!fin.eof()){
       fin>>p0>>p1>>p2>>p3>>p4>>like;
       getline(fin,str);
       p1v.push_back(p1);
       p2v.push_back(p2);
       likev.push_back(like);

       if(p1<p1min) p1min=p1;
       if(p2<p2min) p2min=p2;
       if(p1>p1max) p1max=p1;
       if(p2>p2max) p2max=p2;
       likeTot+=like;
       if(like>likeMax){
	 likeMax=like;
	 p1Best=p1;
	 p2Best=p2;
	 p3Best=p3;
	 p4Best=p4;
       }

    }
    fin.close();

    nentry=(int)(npoint/ndata);
    nz=nentry-6;

    dp2=p2v[1]-p2v[0]; 
    np2=ceil((p2max-p2min)/dp2)+1;
    np1=ndata/np2;
    dp1=(p1max-p1min)/(np1-1);

    cout<<npoint<<"\t"<<ndata<<"\t"<<npoint/ndata<<"\t"<<np1<<"\t"<<np2<<endl;
    cout<<p1min<<"\t"<<p1max<<"\t"<<p2min<<"\t"<<p2max<<"\t"<<dp1<<"\t"<<dp2<<endl;

    //calculate p1 pdf
    cout<<"p1 pdf"<<endl;
    nbins=np1;
    pdf.clear();
    pdf.reserve(nbins);
    cout<<pdf.size()<<endl;
    //for(i=0;i<nbins;i++) pdf[i]=0.0;
    for(i=0;i<nbins;i++) pdf.push_back(0.0);

    xmin=p1min;
    xmax=p1max;
    dx=dp1;

    for(i=0;i<ndata;i++){
       ibin=floor(i/np2);
       pdf[ibin]+=dp2*likev[i]/likeTot;
       //cout<<pdf[ibin]<<endl;
    }

    cout<<pdf.size()<<endl;

    wnorm=0.0;
    for(i=0;i<nbins;i++){
      wnorm+=pdf[i]*dp1;
      //cout<<i<<"\t"<<pdf[i]<<endl;
    }
    //cout<<dp1<<"\t"<<wnorm<<endl;

    files=data_dir+"p1_pdf.dat";
    fout.open(files.c_str());
    for(i=0;i<nbins;i++){
      //cout<<xmin+i*dx<<"\t"<<pdf[i]/wnorm<<"\t"<<wnorm<<"\t"<<likeTot<<endl;
       fout<<xmin+i*dx<<"\t"<<pdf[i]/wnorm<<endl;
    }
    fout.close();

    //calculate p2 pdf
    cout<<"p2 pdf"<<endl;
    nbins=np2;
    pdf.clear();
    pdf.reserve(nbins);
    for(i=0;i<nbins;i++) pdf[i]=0.0;
    xmin=p2min;
    xmax=p2max;
    dx=dp2;

    for(i=0;i<ndata;i++){
       ibin=i%np1;
       pdf[ibin]+=dp1*likev[i]/likeTot;
    }

    wnorm=0.0;
    for(i=0;i<nbins;i++) wnorm+=pdf[i]*dp2;
    
    files=data_dir+"p2_pdf.dat";
    fout.open(files.c_str());
    for(i=0;i<nbins;i++){
       fout<<xmin+i*dx<<"\t"<<pdf[i]/wnorm<<endl;
    }
    fout.close();

       z=8.0;   
    while(z<14.0){
       cout<<z<<endl;
       //calculate xi pdf at given redshift

       iz=floor((z-zmin)/dz)+6;

       files=data_dir+file_like;
       xv.clear();
       fin.open(files.c_str()); 
       while(!fin.eof()){
          for(i=1;i<iz;i++) fin>>trash;
          fin>>xi;
          xv.push_back(xi);
          getline(fin,str);
       }
       fin.close();

       nbins=100;
       pdf.clear();
       pdf.reserve(nbins);
       for(i=0;i<nbins;i++) pdf[i]=0.0;
       xmin=0.0;
       xmax=1.0;
       dx=(xmax-xmin)/(nbins-1);

       for(i=0;i<ndata;i++){
          ibin=floor((xv[i]-xmin)/dx);
          pdf[ibin]+=likev[i]/likeTot/dx;
       }
    
       files=data_dir+"xi_pdf_z"+numToString(z)+".dat";
       fout.open(files.c_str());
       for(i=0;i<nbins;i++){
          fout<<xmin+i*dx<<"\t"<<pdf[i]<<endl;
       }
       fout.close();

       z+=1.0;
    }

    //calculate redshift distribution of xi    
       xi=0.1;
       while(xi<0.99){
          cout<<xi<<endl;

          nbins=40;
          pdf.clear();
          pdf.reserve(nbins);
          for(i=0;i<nbins;i++) pdf[i]=0.0;
          xmin=6.0;
          xmax=20.0;
          dx=(xmax-xmin)/(nbins-1);

          iz=floor((z-zmin)/dz)+6;

          files=data_dir+file_like;
          fin.open(files.c_str()); 
          for(i=0;i<ndata;i++){
             xinew=1.1;
             for(j=1;j<=6;j++) fin>>trash;
             for(j=0;j<nz;j++){
                xiold=xinew;
                fin>>xinew;
                //linearly interpolate between bins
                if(xiold>xi && xi>xinew) z=zmin+(j-1)*dz+dz*(xi-xiold)/(xinew-xiold);
             }
             getline(fin,str);
             ibin=floor((z-xmin)/dx);
             pdf[ibin]+=likev[i]/likeTot/dx;       
          }
          fin.close();

          wnorm=0.0;
          for(i=0;i<nbins;i++) wnorm+=pdf[i]*dx;
          //cout<<xi<<"\t"<<wnorm<<endl;

          files=data_dir+"z_pdf_xi"+numToString(xi)+".dat";
          fout.open(files.c_str());
          for(i=0;i<nbins;i++){
             fout<<xmin+i*dx<<"\t"<<pdf[i]<<endl;
          }
          fout.close();

          xi+=0.1;
       }

       /////////////////////////////////////////////////
       //calculate bestfit history
       /////////////////////////////////////////////////
       double zstep;
       double Nion,gamma,zeta,clump,mfp;


       param.clear();
       param.push_back(p0);
       param.push_back(p1Best);
       param.push_back(p2Best);
       param.push_back(p3Best);
       param.push_back(p4Best);
       Ion.setZetaParam(0.0,param,1);

       files=data_dir+"best_history.dat";
       fout.open(files.c_str());

       fout<<"#"<<p0<<"\t"<<p1Best<<"\t"<<p2Best<<"\t"<<p3Best<<"\t"<<p4Best<<endl;
       fout<<"# like_best like_tot best_prob"<<likeMax<<"\t"<<likeTot<<"\t"<<likeMax/likeTot<<endl;
       fout<<"# tau="<<Ion.getTauCMB()<<endl;
       fout<<"# z Nion  xi  gamma  zeta  mfp  clump"<<endl;
       
       z=20.0;
       zstep=0.5;
       while(z>2.5){
	 Nion=Ion.getNion(z)*MPC*MPC*MPC;
	 xi=Ion.getXI(z);
	 gamma=lyf.getGammaFromNdot(z,Nion);
	 zeta=Ion.getZeta(z);
	 mfp=lyf.mfpFromGamma(z,gamma);
	 clump=Ion.getClumping(z,xi);

	 cout<<z<<"\t"<<Nion<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<mfp<<"\t"<<clump<<endl;
	 fout<<z<<"\t"<<Nion<<"\t"<<xi<<"\t"<<gamma<<"\t"<<zeta<<"\t"<<mfp<<"\t"<<clump<<endl;	
	 z-=zstep;
       }

       fout.close();

///////////////////////////////////////////////////////////////////////
  return 0;
}




