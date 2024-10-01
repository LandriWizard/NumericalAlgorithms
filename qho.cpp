#define _USE_MATH_DEFINES 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "my_header.h"
using namespace std;

#define FALSE 0
#define TRUE 1
#define DEBUG FALSE
#define SESSION 1

#define ARRAYMAX 10000

static double g_E = .5;

void RK4Step(double, double*, void(*) (double, double*, double*), double, int);
void dYdt(double, double*, double*);

double Residual(double);

int main(){

  cout << setiosflags(ios::scientific);
  cout << setprecision(16);
  
  #if DEBUG == TRUE
  cout << "!! DEBUG TIME !!" << endl;
  #endif

  #if SESSION == 1
  double a = 0.;
  double b = 5.;
  int nInters = 11;
  double xL[nInters], xR[nInters];
  double root;
  int nIterations;
  int nRoots;
  Bracket(Residual,a,b,xL,xR,nInters,nRoots);
  int j = 0;
  double tol = 1.e-8;
  double RootsArray[nRoots];
  double code;
  for(int i = 0; i < nRoots; i++){
    code = Bisection(Residual,xL[i],xR[i],tol,tol,nIterations,root);
    if(code == 0){
      RootsArray[j] = root;
      cout << root << endl;
      j++;
    }
  }
  ofstream fdata;
  int N = 800;
  double dE = fabs(a - b)/(double)N;
  double E = 0;
  fdata.open("qho.dat");
  for(int j = 0; j < N; j++){
    fdata << E << " " << Residual(E) << endl;  
    E += dE;
  }
  fdata.close();

  #endif


  return 0;
}

void RK4Step(double t, double* Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){
  double k1[neq], k2[neq], k3[neq], k4[neq];
  double Y1[neq];

  RHSFunc(t,Y,k1);
  for(int i = 0; i < neq; i++){
    Y1[i] = Y[i] + 0.5*dt*k1[i];
  }
  RHSFunc(t+0.5*dt,Y1,k2);
  for(int i = 0; i < neq; i++){
    Y1[i] = Y[i] + 0.5*dt*k2[i];
  }
  RHSFunc(t+0.5*dt,Y1,k3);
  for(int i = 0; i < neq; i++){
    Y1[i] = Y[i] + dt*k3[i];
  }
  RHSFunc(t+dt,Y1,k4);
  for(int i = 0; i < neq; i++){
    Y[i] += 1./6.*dt*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i]);
  }
}

void dYdt(double t, double* Y, double* R){
  double x = Y[0];
  double v = Y[1];
  R[0] = v;
  R[1] = -2.*(g_E - .5*t*t)*x;
}

double Residual(double E){
  g_E = E;
  double xi = -10.;
  double xf = 10.;
  double xm = .8;
  double xl = xi; //clock from left
  double xr = xf; //clock from right
  int N = 800; //points number  
  double w = fabs(xf - xi); //interval width
  double wl = fabs(xm - xi); //left intervals width
  double wr = fabs(xf - xm); //right intervals width
  int NL = int(wl/w*N);
  int NR = int(wr/w*N);
  double dxl = wl/(double)NL;
  double dxr = -wr/(double)NR;
  int neq = 2;
  double YL[neq];
  double YR[neq];
  YL[0] = exp(-.5*xl*xl);
  YL[1] = -xl*exp(-.5*xl*xl);
  YR[0] = exp(-.5*xr*xr);
  YR[1] = -xr*exp(-.5*xr*xr);
  for(int j = 0; j < NL; j++){
    RK4Step(xl,YL,dYdt,dxl,neq);
    xl += dxl;    
  }
  for(int j = 0; j < NR; j++){
    RK4Step(xr,YR,dYdt,dxr,neq);
    xr += dxr;    
  }
  double a = YL[1]*YR[0];
  double b = YR[1]*YL[0];

  return (a - b)/sqrt(a*a + b*b);
}
