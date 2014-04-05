/* multipole.cc
 *
 * spline a data set with x-values x[] and y-values y[] of length n
 *
 */


#include <math.h>
#include <iostream>
#include <fstream>
#include "dnumrecipes.h"
#include "spline2D.h"

using namespace std;


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

Spline2D::Spline2D()
{
  init_flag=0;
  //  cout << "Spline Constructor called successfully." <<endl; 
}

Spline2D::~Spline2D()
{
  //cout <<"Spline destructor has been called." <<endl;

  //if spline initialised then free up vectors
  cleanSpline();

    //cout<<"Spline destructor successful"<<endl;
}

/////////////////////////////////////////////////////////////////////
// Member functions
/////////////////////////////////////////////////////////////////////
int Spline2D::getN()
{
  return n;
}

int Spline2D::getM()
{
  return m;
}

double Spline2D::getX1Max()
{
  return x1max;
}

double Spline2D::getX1Min()
{
  return x1min;
}


double Spline2D::getX2Max()
{
  return x2max;
}

double Spline2D::getX2Min()
{
  return x2min;
}

double Spline2D::getX1element(int i)
{
	if(init_flag==0){
		 cout<<"spline not initiated"<<endl;
		return -1.0;
	}else if(i>m){
		 cout<<"index exceeds spine length"<<endl;
		return -1.0;
	}

	return x1as[i];
}

double Spline2D::getX2element(int i)
{
	if(init_flag==0){
		 cout<<"spline not initiated"<<endl;
		return -1.0;
	}else if(i>n){
		 cout<<"index exceeds spine length"<<endl;
		return -1.0;
	}

	return x2as[i];
}

//return initialisation flag
int Spline2D::checkInit()
{
  return init_flag;
}
/////////////////////////////////////////////////////////////////////
// Utility functions
/////////////////////////////////////////////////////////////////////
void Spline2D::setSplineSP(int m1, int n1, double x1[], double x2[], double **y)
{
  int i, j;

  cleanSpline();

  init_flag=1;
  n=n1;
  m=m1;
  x1as=dvector(1,m);
  x2as=dvector(1,n);
  yas=dmatrix(1,m,1,n);
  y2as=dmatrix(1,m,1,n);  

  //store data
  for(i=1;i<=m;i++) x1as[i]=x1[i];
  for(i=1;i<=n;i++) x2as[i]=x2[i];
  for(i=1;i<=m;i++){
	for(j=1;j<=n;j++){
		yas[i][j]=y[i][j];
	}
  }

  //splineSP(xas,yas,n,natural,natural,y2as);
  splie2(x1as,x2as,yas,m,n,y2as); 

  x1min=x1as[1];
  x1max=x1as[m];
  x2min=x2as[1];
  x2max=x2as[n];

}

// Construct matrix of second derivatives of rows
void Spline2D::splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a)
{
	Spline SP;
	int j;

	for (j=1;j<=m;j++){
		SP.splineSP(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
	}
}


// Construct and evaluate spline to find desired value
void Spline2D::splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y)
{
	int j;
	double *ytmp,*yytmp;
	Spline SP;


	ytmp=dvector(1,m);
	yytmp=dvector(1,m);
	for (j=1;j<=m;j++){
		SP.splintSP(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
	}

	SP.splineSP(x1a,yytmp,m,1.0e30,1.0e30,ytmp);

	SP.splintSP(x1a,yytmp,ytmp,m,x1,y);
	free_dvector(yytmp,1,m);
	free_dvector(ytmp,1,m);
}





/* Constructs spline table for 1-D interpolation.  x is dependent variable
 * vector, y is independent variable vector, and table is returned
 * as y2.  n is the array size.  yp1, yp2 are set to 1.0e30 for natural
 * splines. */
void Spline2D::splineSP(double x[], double y[], int n, double yp1, double ypn,
	    double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;

  u = dvector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0.0;
  else {
    y2[1] = -0.5;
    u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k] = y2[k]*y2[k+1]+u[k];

  free_dvector(u,1,n-1);
  return;
}

/* Spline evaluation function.  xa is the dependent variable vector
 * and ya is the independent variable vector.  y2a is the spline table
 * from a previous call to spline.  n is the array size.  x is the point
 * at which the interpolation is to be performed.  Result is returned
 * in y. NOTE xa must be in increasing order! */
void Spline2D::splintSP(double xa[], double ya[], double y2a[], int n, double x,
	    double *y)
{
  int klo,khi,k;
  double h,b,a;
  static int kloOld,khiOld,useflag;

  if(useflag!=0 && (x<xa[khiOld] && x>xa[khiOld])){
    khi=khiOld;
    klo=kloOld;
    //  }else if(useflag!=0){
    //hunt(xa,n,x,&klo);
    //khi=klo+1;
  }else{
    //resort to bisection
    klo=1;
    khi=n;

    while (khi-klo > 1) {
      k = (khi + klo) >> 1;
      if (xa[k] > x)
	khi = k;
      else
	klo = k;
    }
  }
  // khi and klo now bracket the input value of x
  h = xa[khi] - xa[klo];
  if (h == 0.0)
    throw NumRecException("Bad xa input into routine splint");
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;
  *y = a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + 
				(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  kloOld=klo;
  khiOld=khi;
  useflag=1;
  return;
}

// Dummy function to return splined value
double Spline2D::returnValue(double x1, double x2)
{
  double y;
  if(x1<x1min) {
    cout<<"x1 value below range of spline: "<<x1<<"\t"<<x1min<<"\t"<<x2<<"\t"<<x2min<<endl;
    return 0.0;
  } else if(x1>x1max){
    cout<<"x1 value above range of spline: "<<x1<<"\t"<<x1max<<"\t"<<x2<<"\t"<<x2max<<endl;
    return 0.0;
  }

  if(x2<x2min) {
    cout<<"x2 value below range of spline: "<<x1<<"\t"<<x1min<<"\t"<<x2<<"\t"<<x2min<<endl;
    return 0.0;
  } else if(x2>x2max){
    cout<<"x2 value above range of spline: "<<x2<<"\t"<<x2max<<endl;
    return 0.0;
  }


  splin2(x1as,x2as,yas,y2as,m,n,x1,x2,&y);
  return y;
}


//bisection is not the best method for large vectors
//instead use hunt to find the desired limits
void Spline2D::hunt(double xx[], int n, double x, int *jlo)
{
	int jm,jhi,inc;
	int ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo)--;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc <<= 1;
				if (inc >= jhi) {
					*jlo=0;
					break;
				}
				else *jlo=jhi-inc;
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}

//////////////////////////////////////////////////////////////////////
//  Load spline from file
//////////////////////////////////////////////////////////////////////
//

// Load spline information from a file with
// first line giving m and n
// then m*n lines giving (x1,x2,y) triples
//
void Spline2D::loadFileSpline(string file)
{
	double trash, *x1, *x2, **y;
	int i,j;
	int m, n, nrow;
	ifstream fin;
	char buffer[100];

	fin.open(file.c_str());
	if(fin.fail()){
		 cout<<"File not found"<<endl;
		 return;
	}

	// read in m and n values
	fin>>m>>n;

	nrow=0;
	while(!fin.eof()){
		fin.getline(buffer,100);
		nrow++;
	}
	nrow--; //loop goes one place too far
	nrow--; //ignore line containing m,n
	fin.close();


	//cout<<m<<"\t"<<n<<"\t"<<m*n<<"\t"<<nrow<<endl;

	//check we have correct number of data points
	if(m*n!=nrow){
		cout<<"error dimensions and number of data entries do not agree"<<endl;
	return;
	}

	x1=dvector(1,m);
	x2=dvector(1,n);
	y=dmatrix(1,m,1,n);
	fin.open(file.c_str());

	fin>>m>>n;

	for(i=1;i<=m;i++){
		for(j=1;j<=n;j++){
			fin>>x1[i]>>x2[j]>>y[i][j];
			//cout<<x1[i]<<"\t"<<x2[j]<<"\t"<<y[i][j]<<endl;
		}
	}
	fin.close();
			
	setSplineSP(m,n,x1,x2,y);

	free_dvector(x1,1,m);
	free_dvector(x2,1,n);
	free_dmatrix(y,1,m,1,n);
}


///////////////////////////////////////////////////////////////
// Clean up spline for reuse
///////////////////////////////////////////////////////////////
void Spline2D::cleanSpline()
{

    if(init_flag==1){
      free_dvector(x1as,1,m);
      free_dvector(x2as,1,n);
      free_dmatrix(yas,1,m,1,n);
      free_dmatrix(y2as,1,m,1,n);
    }

    init_flag=0;
    n=0;
    m=0;
    x1min=0.0;
    x1max=0.0;
    x2min=0.0;
    x2max=0.0;
}

////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////
