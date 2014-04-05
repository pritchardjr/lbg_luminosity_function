
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
     int like_calc_flag(1);
     int nthread(0);
     int gamma_flag(0);
     int sort_flag(0);
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
    case 's':
      sort_flag=atoi(&argv[1][2]);
      break;
    case 'p':
      param_flag=atoi(&argv[1][2]);
      break;
    case 'c':
     cmb_flag=atoi(&argv[1][2]);
      break;
    case 'm':
     tocm_flag=atoi(&argv[1][2]);
      break;
    case 'n':
     ndot_flag=atoi(&argv[1][2]);
      break;

    default:
      cerr << "Bad option" <<argv[1] <<"\n";
    }
    --argc;
    ++argv;
  }

///////////////////////////////////////////////////////////////////////
    // Post process cumulative frequency tables to calculate contour plot
///////////////////////////////////////////////////////////////////////
    string files;
    int nz(0), nx(0), ndata(0);
    double z, x, cpdf, trash;
    double *zv, *xv, *cpdfv;
    double zmax(-1e30),xmax(-1e30);
    Spline pdfSP;
    int il;
    string dirbase;
    string prefix;

    dirbase="./ndot/";
    if(param_flag==2)     dirbase="./twostep/";

    if(cmb_flag==0 && ndot_flag==0) dirbase=dirbase+"bolton/";
    if(cmb_flag==1 && ndot_flag==0) dirbase=dirbase+"wmap3/";
    if(cmb_flag==2 && ndot_flag==0) dirbase=dirbase+"planck/";

    if(ndot_flag==2) dirbase=dirbase+"nolya/";
    if(ndot_flag==3) dirbase=dirbase+"betterg/";

    if(tocm_flag==1)  dirbase=dirbase+"21cm/";

    if(sort_flag==0) prefix="xi";
    if(sort_flag==1) prefix="sn";
    if(sort_flag==2) prefix="dxdz";
    files=dirbase+prefix+"_contour.dat";
    cout<<files<<endl;
    fin.open(files.c_str());
    while(!fin.eof()){
      if(sort_flag==0){
	fin>>z>>x>>cpdf>>trash>>trash;
      }else{
	fin>>z>>x>>cpdf>>trash;	
      }
      //cout<<z<<"\t"<<x<<"\t"<<cpdf<<endl;
      ndata++;
      if(x>xmax){
	xmax=x;
	nx++;
      }
      if(z>zmax){
	zmax=z;
	nz++;
      }
    }
    ndata--;  //goes one space too far
    fin.close();

    cout<<nz<<"\t"<<nx<<"\t"<<nz*nx<<"\t"<<ndata<<endl;

    zv=dvector(1,nx);
    xv=dvector(1,nx);
    cpdfv=dvector(1,nx);

    //one and two sigma levels
    vector<double> lev;
    lev.push_back(0.0228);
    lev.push_back(0.159);
    lev.push_back(0.5);
    lev.push_back(0.841);
    lev.push_back(0.9773);

    

    //now by redshift create spline and then find the levels
    fin.open(files.c_str());
    files=dirbase+prefix+"_contourprep.dat";
    fout.open(files.c_str());

    //redshift steps
    for(i=1;i<=nz;i++){
      //prepare spline for this redshift
      for(j=1;j<=nx;j++){
	if(sort_flag==0) fin>>zv[j]>>xv[j]>>cpdfv[j]>>trash>>trash;
	if(sort_flag==1) fin>>zv[j]>>trash>>cpdfv[j]>>xv[j];	
	if(sort_flag==2) fin>>zv[j]>>xv[j]>>cpdfv[j]>>trash;
	}
      //for(j=1;j<=nx;j++) cout<<zv[j]<<"\t"<<xv[j]<<"\t"<<cpdfv[j]<<endl;      
      
      pdfSP.setSplineSP(nx,xv,cpdfv);
      //cout<<"spline good"<<endl;

      fout<<zv[1];
      cout<<zv[1]<<endl;
      //use spline to find x corresponding to lev
      for(il=0;il<=lev.size()-1;il++){
	//handle case where lev below lower limit of binning
	if(lev[il]<pdfSP.returnValue(pdfSP.getXMin())){
	  x=pdfSP.getXMin();
	}else{
	  //cout<<lev[il]<<"\t"<<lev.size()<<endl;
	  x=pdfSP.returnOrdinate(lev[il]);
	}
	fout<<"\t"<<x;
      }
      fout<<endl;

    }

///////////////////////////////////////////////////////////////////////
  return 0;
}




