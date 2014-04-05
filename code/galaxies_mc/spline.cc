/* multipole.cc
 *
 * spline a data set with x-values x[] and y-values y[] of length n
 *
 */


#include <math.h>
#include <iostream>
#include <fstream>
#include "dnumrecipes.h"
#include "spline.h"

using namespace std;


/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

Spline::Spline()
{
  init_flag=0;
  //  cout << "Spline Constructor called successfully." <<endl; 
}

Spline::~Spline()
{
  //cout <<"Spline destructor has been called." <<endl;

  //if spline initialised then free up vectors
  cleanSpline();

    //cout<<"Spline destructor successful"<<endl;
}

/////////////////////////////////////////////////////////////////////
// Member functions
/////////////////////////////////////////////////////////////////////
int Spline::getN()
{
  return n;
}

double Spline::getXMax()
{
  return xmax;
}

double Spline::getXMin()
{
  return xmin;
}

double Spline::getXelement(int i)
{
	if(init_flag==0){
		 cout<<"spline not initiated"<<endl;
		return -1.0;
	}else if(i>n){
		 cout<<"index exceeds spine length"<<endl;
		return -1.0;
	}

	return xas[i];
}

//return initialisation flag
int Spline::checkInit()
{
  return init_flag;
}
/////////////////////////////////////////////////////////////////////
// Utility functions
/////////////////////////////////////////////////////////////////////
void Spline::setSplineSP(int n1, double x[], double y[])
{
  int i;
  double natural(1.0e30);

  //first clean up old spline if necessary
  cleanSpline();

  init_flag=1;
  n=n1;
  xas=dvector(1,n);
  yas=dvector(1,n);
  y2as=dvector(1,n);  

  for(i=1;i<=n;i++){
    xas[i]=x[i];
    yas[i]=y[i];
   }
  splineSP(xas,yas,n,natural,natural,y2as);
 
  //   cout<<xas[1]<<"\t"<<yas[1]<<"\t"<<y2as[1]<<endl;
  //for(i=2;i<=n;i++){
    // cout<<xas[i]<<"\t"<<yas[i]<<"\t"<<y2as[i]<<"\t"<<xas[i]-xas[i-1]<<endl;
  //}

  xmin=xas[1];
  xmax=xas[n];
}


/* Constructs spline table for 1-D interpolation.  x is dependent variable
 * vector, y is independent variable vector, and table is returned
 * as y2.  n is the array size.  yp1, yp2 are set to 1.0e30 for natural
 * splines. */
void Spline::splineSP(double x[], double y[], int n, double yp1, double ypn,
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
void Spline::splintSP(double xa[], double ya[], double y2a[], int n, double x,
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
double Spline::returnValue(double x)
{
  double y;
  if(x<xmin) {
    cout<<"x value below range of spline: "<<x<<"\t"<<xmin<<endl;
    return 0.0;
  } else if(x>xmax){
    cout<<"x value above range of spline: "<<x<<"\t"<<xmax<<endl;
    return 0.0;
  }


  splintSP(xas,yas,y2as,n,x,&y);
  return y;
}


//bisection is not the best method for large vectors
//instead use hunt to find the desired limits
void Spline::hunt(double xx[], int n, double x, int *jlo)
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
// Load spline information from a file with ncol columns which has
// the desired x information in the xcol th column and y info in
// the ycol th column
//
void Spline::loadFileSpline(string file, int ncol, int xcol, int ycol)
{
	double trash, *x, *y;
	int i,j,nrow;
	ifstream fin;
	char buffer[100];

//	cout<<file<<endl;

	fin.open(file.c_str());
	if(fin.fail()){
		 cout<<"File not found"<<endl;
		 return;
	}

	nrow=0;
	while(!fin.eof()){
		fin.getline(buffer,100);
		nrow++;
	}
	nrow--; //loop goes one place too far
	fin.close();
//	cout<<nrow<<endl;

	x=dvector(1,nrow);
	y=dvector(1,nrow);
	fin.open(file.c_str());
	for(i=1;i<=nrow;i++){
		for(j=1;j<=ncol;j++){
			if(j==xcol){
				 fin>>x[i];
			}else if(j==ycol){
				 fin>>y[i];
			}else {
				fin>>trash;
			}
		}
//		cout<<x[i]<<"\t"<<y[i]<<endl;
	}
	fin.close();
			
	setSplineSP(nrow,x,y);

	free_dvector(x,1,nrow);
	free_dvector(y,1,nrow);
}

///////////////////////////////////////////////////////////////
// Clean up spline for reuse
///////////////////////////////////////////////////////////////
void Spline::cleanSpline()
{

    if(init_flag==1){
      free_dvector(xas,1,n);
      free_dvector(yas,1,n);
      free_dvector(y2as,1,n);
    }

    init_flag=0;
    n=0;
    xmin=0.0;
    xmax=0.0;
}

///////////////////////////////////////////////////////////////
// spline inversion
///////////////////////////////////////////////////////////////

//use interpolation to get x value corresponding to y
double Spline::returnOrdinate(double y)
{
  double ymin(1.0e30), ymax(-1.0e30), yuse,xuse;
  double tol(1.0e-5);
  int i;

  //check y lies within the range of the spline
  //THIS SHOULD GO IN SETSPLINESP
  for(i=1;i<=n;i++){
    yuse=yas[i];
    ymin=fmin(yuse,ymin);
    ymax=fmax(yuse,ymax);
  }
  //cout<<ymin<<"\t"<<ymax<<endl;
  if(y<ymin){
    cout<<"y value below spline range: "<<y<<"\t"<<ymin<<endl;
    return getXMin();
  }
  if(y>ymax){
    cout<<"y value above spline range: "<<y<<"\t"<<ymax<<endl;
    return getXMax();
  }

  //now find the value of x corresponding to y
  //cout<<"root finding"<<xmin<<"\t"<<xmax<<"\t"<<y<<endl;
  xuse=zriddrConst(y,xmin,xmax,tol);

  return xuse;
}

// Same as zriddr, but function isn't expected to calculate its derivative
double Spline::zriddrConst(double yval, double x1, double x2, double xacc)
{ 
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

  fl = returnValue(x1) - yval;
  fh = returnValue(x2) - yval;
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl = x1;
    xh = x2;
    ans = -1.11e30;
    for (j=1;j<=MAXIT;j++) {
      xm = 0.5*(xl + xh);
      fm = returnValue(xm) - yval;
      s = sqrt(fm*fm - fl*fh);
      if (s == 0.0)
	return ans;
      xnew = xm + (xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xacc) 
	return ans;
      ans = xnew;
      fnew = returnValue(ans) - yval;
      if (fnew == 0.0) 
	return ans;
      if (sign(fm,fnew) != fm) {
	xl = xm;
	fl = fm;
	xh = ans;
	fh = fnew;
      } else {
	if (sign(fl,fnew) != fl) {
	  xh = ans;
	  fh = fnew;
	} else {
	  if (sign(fh,fnew) != fh) {
	    xl = ans;
	    fl = fnew;
	  } else {
	    throw NumRecException("never get here");
	  }
	}
      }
      if (fabs(xh-xl) <= xacc)
	return ans;
    }
    throw NumRecException("Maximum number of iterations exceeded in zriddr");
  } else {
    if (fl == 0.0) 
      return x1;
    if (fh == 0.0)
      return x2;
    throw NumRecException("root must be bracketed in zriddr");
  }
  return 0.0;
}

