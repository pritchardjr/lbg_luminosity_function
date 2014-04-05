/* atomic.cc
 * Contains functions for calculating useful spectroscopic quantites
 *
  */

#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include "atomic.h"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_result.h>
#include "dnumrecipes.h"


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

Atomics::Atomics(int rNMAX_in, double **pastore, double **pastore_thick)
{
  // cout <<"Atomic constructor has been called." <<endl;
  
  thick_flag=0;
  rNMAX=rNMAX_in;
  //pastore=dmatrix(1,rNMAX,1,2*(rNMAX-1)-1);

  initRecycleLyman(pastore);
  
  pstore=dmatrix(1,rNMAX,1,2*(rNMAX-1)-1);

   //   cout <<pstore[2][3]<<"\t"<<pastore[2][3]<<endl;
  //pastore_thick=dmatrix(1,rNMAX,1,2*(rNMAX-1)-1);
  initRecycleLymanThick(pastore_thick);


  //cout << "Atomic Constructor called successfully." <<endl; 
  
}

Atomics::~Atomics()
{
  // cout <<"Atomic destructor has been called." <<endl;
  //free_dmatrix(pastore,1,rNMAX,1,2*(rNMAX-1)-1);
  //free_dmatrix(pastore_thick,1,rNMAX,1,2*(rNMAX-1)-1);
}

/*********************************************************************
 **************** Utility Functions      *****************************
 ********************************************************************/
// Calculates the transition probability A_{ul} between an upper 
// state of hydrogen  u=(n,l,j) and lower state l=(n',l',j').  
// 
double Atomics::transitionA(double n, double l, double j, double np, double lp,
			    double jp)
{
  double tranA(1.0);
  double lmax;
  double coupling;

  coupling=gsl_sf_coupling_6j(int(2.0*l),int(2.0*j),int(2.0*0.5),int(2.0*jp),int(2.0*lp),int(2.0*1.0));

  // In optically thick case transitions to 1S state are regenerated
  if(thick_flag ==1 && ((np==1.0)&&(lp==0.0))) return 0.0;

  if(coupling == 0.0){
    return 0.0;
  }

  lmax=fmax(fabs(l),fabs(lp));
  tranA = lmax*(2.0*jp+1.0)*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS;
  tranA *= BOHRA0*BOHRA0*coupling*coupling;
  tranA *= pow(hydrogenMatrix(n,l,np,lp),2.0);
  tranA *= 64.0*PI*PI*PI*PI/3.0/PLANCKCONSTANT;
  tranA /= pow(hydrogenWave(n,np),3.0);
  
  return tranA;
}

//Need to modify below so that it works for (n,l) to (np,lp) transitions.

// Calculates the absorption oscillator strength of the transition  between
// an upper state of hydrogen  u=(n,l,j) and lower state l=(n',l',j')
double Atomics::transitionF(double n, double l, double j, double np, double lp,
			    double jp)
{
  double transF;

  transF = ELECTRONMASS*pow(SPEEDOFLIGHT_CGS,3.0);
  transF/= 8.0*PI*PI*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS*pow(hydrogenFreq(n,np),2.0);
  transF *= transitionA(n,l,j,np,lp,jp);
  transF *= (2.0*j+1.0)/(2.0*jp+1.0);

  return transF;
} 

// Calculates the frequency of the transition between the upper level n
// and lower level np in the hydrogen atom.
double Atomics::hydrogenFreq(double n, double np)
{
  double hf(1.0);
  
  hf=RYDBERG_CGS*(1.0/np/np - 1.0/n/n)/PLANCKCONSTANT;

  return hf;
}

// Calculates the wavelength of the transition between the upper level n
// and lower level np in the hydrogen atom.
double Atomics::hydrogenWave(double n, double np)
{
  double hl(SPEEDOFLIGHT_CGS);
  
  hl /= RYDBERG_CGS*(1.0/np/np - 1.0/n/n)/PLANCKCONSTANT;
  
  return hl;
}

// Calculates the Energy of the transition between the upper level n
// and lower level np in the hydrogen atom.
double Atomics::hydrogenEnergy(double n, double np)
{
  double he(RYDBERG_CGS);
  
  he *= (1.0/np/np - 1.0/n/n);
  
  return he;
}


// Calculates the transition strength S_{ul} between an upper 
// state of hydrogen  u=(n,l,j) and lower state l=(n',l',j').  
//By convention the 
// 
double Atomics::transitionS(double n, double l, double j, double np, double lp,
			    double jp)
{
  double tranS(1.0);
  double lmax;
  double coupling;

  coupling=gsl_sf_coupling_6j(int(2.0*l),int(2.0*j),int(2.0*0.5),int(2.0*jp),int(2.0*lp),int(2.0*1.0));

  if(coupling == 0.0){
    return 0.0;
  }

  lmax=fmax(fabs(l),fabs(lp));
  tranS = lmax*(2.0*j+1.0)*(2.0*jp+1.0)*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS;
  tranS *= BOHRA0*BOHRA0*coupling*coupling;
  tranS *=pow(hydrogenMatrix(n,l,np,lp),2.0);
  
  return tranS;
}

// Calculates the total decay rate from a state u=(n,l,j).
//
// Loop structure could probably be made more elegant
double Atomics::transitionATotal(double n, double l, double j)
{
  double Atotal(0.0);
  double m,jp;
  
  // For l>0 case handle transitions which have Delta l = -1
  if(l > 0.0){
    m=l;
    while(m<=n-1.0){
      jp=fabs(l-1.5);
      while(jp<=l-0.5){
	Atotal+=transitionA(n,l,j,m,l-1.0,jp);
	jp += 1.0;
      }
      	m += 1.0;
    }
  }
  // For all cases handle transitions with Delta l = +1
  m=l+2.0;
  while(m<=n-1.0){
    jp=fabs(l+0.5);
    while(jp<=l+1.5){
      Atotal+=transitionA(n,l,j,m,l+1.0,jp);
      jp += 1.0;
    }
    m += 1.0;
  } 
 
  return Atotal;

}

// Calculates the probability P_ul that a photon in the state u=(n,l,j) will
// make a transition to the state l=(n',l',j')
  double Atomics::transitionP(double n, double l, double j, double np,
			      double lp, double jp)
{
  double transP(1.0);

  transP=transitionA(n,l,j,np,lp,jp)/transitionATotal(n,l,j);

  return transP;
}


// Calculates the matrix element R^{n',l'}_{n,l}/a_0 for the hydrogen 
// atom.  Here (n',l') describe the lower level and (n,l) the upper level.
// a_0 is the Bohr radius.
// Reference: Rudnick (1935) Physical Review 48, 807
double Atomics::hydrogenMatrix(double n, double l, double np, double lp)
{
  double matR;

  matR=factorial(n+l)/factorial(np+lp)/factorial(n-l-1.0)/factorial(np-lp-1.0);
  matR=sqrt(matR);
  matR *= pow(n-np,n-np-1.0)/pow(n+np,n+np+1.0);
  matR *= pow(2.0,l+lp+4.0)*pow(n,lp+3.0)*pow(np,l+3.0);
  matR *= polyR(n,l,np,lp);
  return matR;
}

// Calculates the polynomial required by hydrogenMatrix 
double Atomics::polyR(double n, double l, double np, double lp)
{
  double pR;
  double lmax,v,nr,r;

  lmax=fmax(fabs(l),fabs(lp));
  r=np-lmax-1.0;
  nr=n-l-1.0;
  v=-4.0*n*np/pow(n-np,2.0);

  if(fabs(lp-l+1.0)<1.0e-5){
    //This is the delta l = -1 case
    pR = (n+np)*term2F1(-r,-n+lmax+1.0,2*lmax+1.0,v);
    pR -= (n-np)*term2F1(-r,-n+lmax,2*lmax+1.0,v);
    pR *=pow(-1.0,r)*pow(n-np,2.0*r)*factorial(2.0*lmax+r)/factorial(2.0*lmax);
    pR /= 2.0*np;
    return pR;
  }

   if(fabs(lp-l-1.0)<1.0e-5){
    //This is the delta l = +1 case
    pR = (n+np)*(n-lmax)*term2F1(-r,-n+lmax+1.0,2*lmax+1.0,v);
    pR -= (n-np)*(n+lmax)*term2F1(-r,-n+lmax,2*lmax+1.0,v);
    pR *=pow(-1.0,r)*pow(n-np,2.0*r)*factorial(2.0*lmax+r)/factorial(2.0*lmax);
    pR /= 2.0*n;
    return pR;
  } 

   cout << "Delta l not equal to +- 1 in PolyR!\n" ;
   return 0.0;
}

// Calculate the lyman alpha recycling fractions for the state u=(n,l,j)
//
double Atomics::recycleLyman(double nin, double lin, double jin,
			    int nmax)
{
  double recycleL(0.0);
  double **pa;
  double m,jp;
  double n,l,j;
  int p,q,r;
  
  // if nmax less than 4 no point in running the code
  if(nmax < 4){
    cout << "nmax is less than 4 in recycleLyman\n";
    return 0.0;
  }

  // initialise storage matrix. Each level has 2 spin states +-1/2
  // store by (n, 2*l+(j-l+0.5) = (n,l+j+1/2)
  pa=dmatrix(1,nmax,1,2*(nmax-1)-1);

  // Specify the easily determined recycling fractions
  pa[1][1]=0.0;
  pa[2][1]=0.0;
  pa[2][2]=1.0;
  pa[2][3]=1.0;
  pa[3][1]=1.0;
  pa[3][2]=0.0;
  pa[3][3]=0.0;
  pa[3][4]=1.0;
  pa[3][5]=1.0;

  // Now calculate recycling fractions iteratively from n=4 up to nmax
  // Uses same loop as transitionATotal


  for(p=4;p<=nmax;p++){
    n=double(p);
    for(q=0;q<=p-1;q++){
      l=double(q);
      for(r=0;r<2;r++){
	if(l==0.0){
	  j=0.5;
	  r=1;
	}else{
	  j=l-0.5+double(r);
	} 

	recycleL=0.0;
	// Following code calculates for an individual level

	// For l>0 case handle transitions which have Delta l = -1
	if(l > 0.0){
	  m=l;
	  while(m<=n-1.0){
	    jp=fabs(l-1.5);
	    while(jp<=l-0.5){
	      // Need to store both j states l-1/2, l+1/2
	      recycleL+=transitionP(n,l,j,m,l-1.0,jp)*pa[int(m)][int(l+jp-0.5)];
	      jp += 1.0;
	    }
	    m += 1.0;
	  }
	}
	// For all cases handle transitions with Delta l = +1
	m=l+2.0;
	while(m<=n-1.0){
	  jp=fabs(l+0.5);
	  while(jp<=l+1.5){
	    recycleL+=transitionP(n,l,j,m,l+1.0,jp)*pa[int(m)][int(l+jp+1.5)];
	    jp += 1.0;
	  }
	  m += 1.0;
	} 

	pa[p][2*q+r]=recycleL;
      }
    }
  }
  return pa[int(nin)][int(lin+jin+0.5)];
}

  

// Initiate the matrix pastore with the Lyman alpha recycling fractions
// for levels up to n=rNMAX.  Optically thin case.
double Atomics::initRecycleLyman(double **pastore)
{
  double recycleL(0.0);
  double m,jp;
  double n,l,j;
  int p,q,r;

  int reset_flag(0);
  int rNMAX_file;
  char *file="./pastore_thin.dat";
  ifstream fin(file);
  ofstream fout;

  thick_flag=0;

  // if nmax less than 4 no point in running the code
  if(rNMAX < 4){
    cout << "nmax is less than 4 in recycleLyman\n";
    return 0.0;
  }

  //If it exists and is big enough then use file data
  if(fin){
    fin >> rNMAX_file;
    if(rNMAX<=rNMAX_file){
      // cout <<"reading data"<<endl;
      // read in data
      fin >>pastore[1][1];
      for(p=2;p<=rNMAX;p++){
	for(q=0;q<=p-1;q++){
	  for(r=0;r<2;r++){
	    fin >> pastore[p][2*q+r];
	  }
	}
      }
      fin.close();
      thick_flag=0;
      return 1.0;
    }
    // not enough data in file. Must rest it.
    reset_flag=1;
  }
  //If it doesn't exist then calculate from scratch
  if((!fin) || (reset_flag==1)){
    // use pastore storage matrix. Each level has 2 spin states +-1/2
    // store by (n, 2*l+(j-l+0.5) = (n,l+j+1/2)

    // Specify the easily determined recycling fractions
    pastore[1][1]=0.0;
    pastore[2][1]=0.0;
    pastore[2][2]=1.0;
    pastore[2][3]=1.0;
    pastore[3][1]=1.0;
    pastore[3][2]=0.0;
    pastore[3][3]=0.0;
    pastore[3][4]=1.0;
    pastore[3][5]=1.0;
    
    // Now calculate recycling fractions iteratively from n=4 up to nmax
    // Uses same loop as transitionATotal
    
    
    for(p=4;p<=rNMAX;p++){
      n=double(p);
      for(q=0;q<=p-1;q++){
	l=double(q);
	for(r=0;r<2;r++){
	  if(l==0.0){
	    j=0.5;
	    r=1;
	  }else{
	    j=l-0.5+double(r);
	  } 
	  
	  recycleL=0.0;
	  // Following code calculates for an individual level
	  
	  // For l>0 case handle transitions which have Delta l = -1
	  if(l > 0.0){
	    m=l;
	    while(m<=n-1.0){
	      jp=fabs(l-1.5);
	      while(jp<=l-0.5){
		// Need to store both j states l-1/2, l+1/2
		recycleL+=transitionP(n,l,j,m,l-1.0,jp)*pastore[int(m)][int(l+jp-0.5)];
		jp += 1.0;
	      }
	      m += 1.0;
	    }
	  }
	  // For all cases handle transitions with Delta l = +1
	  m=l+2.0;
	  while(m<=n-1.0){
	    jp=fabs(l+0.5);
	    while(jp<=l+1.5){
	      recycleL+=transitionP(n,l,j,m,l+1.0,jp)*pastore[int(m)][int(l+jp+1.5)];
	      jp += 1.0;
	    }
	    m += 1.0;
	  } 
	  
	  pastore[p][2*q+r]=recycleL;
	}
      }
    }
  
    //store this data in a file
    //cout <<"writing data"<<endl;
    fout.open(file);
    fout << rNMAX<<endl;
    // read in data
    fout <<pastore[1][1]<<endl;
    for(p=2;p<=rNMAX;p++){
      for(q=0;q<=p-1;q++){
	for(r=0;r<2;r++){
	  fout << pastore[p][2*q+r]<<endl;
	}
      }
    }
    fout.close();

    thick_flag=0;
    return 1.0;
  }

  cout << "Error in initLymanRecycleThin"<<endl;
  return 0.0;
}

// Initiate the matrix pastore with the Lyman alpha recycling fractions
// for levels up to n=rNMAX.  Optically thick case.
double Atomics::initRecycleLymanThick(double **pastore_thick)
{
  double recycleL(0.0);
  double m,jp;
  double n,l,j;
  int p,q,r;

  int reset_flag(0);
  int rNMAX_file;
  char *file="./pastore_thick.dat";
  ifstream fin(file);
  ofstream fout;

  thick_flag=1;

  // if nmax less than 4 no point in running the code
  if(rNMAX < 4){
    cout << "nmax is less than 4 in recycleLyman\n";
    return 0.0;
  }

  //If it exists and is big enough then use file data
  if(fin){
    fin >> rNMAX_file;
    if(rNMAX<=rNMAX_file){
      //cout <<"reading data"<<endl;
      // read in data
      fin >>pastore_thick[1][1];
      for(p=2;p<=rNMAX;p++){
	for(q=0;q<=p-1;q++){
	  for(r=0;r<2;r++){
	    fin >> pastore_thick[p][2*q+r];
	  }
	}
      }
      fin.close();
      thick_flag=0;
      return 1.0;
    }
    // not enough data in file. Must rest it.
    reset_flag=1;
  }
  //If it doesn't exist then calculate from scratch
  if((!fin) || (reset_flag==1)){
    // use pastore storage matrix. Each level has 2 spin states +-1/2
    // store by (n, 2*l+(j-l+0.5) = (n,l+j+1/2)

    // Specify the easily determined recycling fractions
    pastore_thick[1][1]=0.0;
    pastore_thick[2][1]=0.0;
    pastore_thick[2][2]=1.0;
    pastore_thick[2][3]=1.0;
    pastore_thick[3][1]=1.0;
    pastore_thick[3][2]=0.0;
    pastore_thick[3][3]=0.0;
    pastore_thick[3][4]=1.0;
    pastore_thick[3][5]=1.0;
    
    // Now calculate recycling fractions iteratively from n=4 up to nmax
    // Uses same loop as transitionATotal
    
    
    for(p=4;p<=rNMAX;p++){
      n=double(p);
      for(q=0;q<=p-1;q++){
	l=double(q);
	for(r=0;r<2;r++){
	  if(l==0.0){
	    j=0.5;
	    r=1;
	  }else{
	    j=l-0.5+double(r);
	  } 
	  
	  recycleL=0.0;
	  // Following code calculates for an individual level
	  
	  // For l>0 case handle transitions which have Delta l = -1
	  if(l > 0.0){
	    m=l;
	    while(m<=n-1.0){
	      jp=fabs(l-1.5);
	      while(jp<=l-0.5){
		// Need to store both j states l-1/2, l+1/2
		recycleL+=transitionP(n,l,j,m,l-1.0,jp)*pastore_thick[int(m)][int(l+jp-0.5)];
		jp += 1.0;
	      }
	      m += 1.0;
	    }
	  }
	  // For all cases handle transitions with Delta l = +1
	  m=l+2.0;
	  while(m<=n-1.0){
	    jp=fabs(l+0.5);
	    while(jp<=l+1.5){
	      recycleL+=transitionP(n,l,j,m,l+1.0,jp)*pastore_thick[int(m)][int(l+jp+1.5)];
	      jp += 1.0;
	    }
	    m += 1.0;
	  } 
	  
	  pastore_thick[p][2*q+r]=recycleL;
	}
      }
    }
  
    //store this data in a file
    //cout <<"writing data"<<endl;
    fout.open(file);
    fout << rNMAX<<endl;
    // read in data
    fout <<pastore_thick[1][1]<<endl;
    for(p=2;p<=rNMAX;p++){
      for(q=0;q<=p-1;q++){
	for(r=0;r<2;r++){
	  fout << pastore_thick[p][2*q+r]<<endl;
	}
      }
    }
    fout.close();

    thick_flag=0;
    return 1.0;
  }

  cout << "Error in initLymanRecycleThick"<<endl;
  return 0.0;
}


/*******************************************************
 ************* Data retrieval functions ******************
 ********************************************************/

  double Atomics::getRecycleLyman(double n, double l, double j, double **pastore)
{

  return pastore[int(n)][int(l+j+0.5)];
}

  double Atomics::getRecycleLymanThick(double n, double l, double j, double **pastore_thick)
{

  return pastore_thick[int(n)][int(l+j+0.5)];
}


/*********************************************************************
 ********* Math functions for Hypergeometric Functions   *************
 ********************************************************************/
//Calculates the value of the Pochhammer symbol (a)_k
double Atomics::intPochhammer(double a, double k)
{
  double temp(1.0);
  double j(0.0);

  if (k <0.0) {
    cout << "k negative in intPochhammer!\n";
    return 0.0;
  }
  if(k ==0.0) {
    return 1.0;
  }

  while(j <= k-1.0){
    temp *= (a+j);
    j +=1.0;
      }

  return temp;
}

//Calculate the factorial n!
double Atomics::factorial(double n)
{
  double temp(1.0);
  double j(1.0);

  if(n<0.0){
    cout << "negative argument to factorial!\n";
    return 0.0;
  }

  while(j<n){
    j += 1.0;
    temp *= j;
  }
  return temp;
}
	

// Calculates the hypergeometric function 2F1(a,b;c,x)
// Requires that a or b be a negative integer so that the series
// terminates
double Atomics::term2F1(double a, double b, double c, double x)
{
  double temp(1.0);
  double k(1.0);
  double kmax(0.0);

  // Check to make sure that at least one of a or b are negative
  if( (a>0.0) && (b>0.0)){
    cout << "a and b both positive in term2F1. Series must terminate\n";
    return 0.0;
  }
  // Work out when series terminates
  if( (a>0.0) && (b<0.0)) {
    kmax=-b;
  }
   if( (a<0.0) && (b>0.0)) {
    kmax=-a;
  } 
   if( (a<0.0) && (b<0.0)) {
    kmax=fmin(fabs(a),fabs(b));
  } 

   while( k<= kmax){
     temp += intPochhammer(a,k)*intPochhammer(b,k)*pow(x,k)/intPochhammer(c,k)/factorial(k);
     k +=1.0;
   }
   return temp;
}






    

