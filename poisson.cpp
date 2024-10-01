#define _USE_MATH_DEFINES 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "my_header.h"
using namespace std;

#define FALSE 0
#define TRUE 1
#define DEBUG TRUE
#define SESSION 1

#define ARRAYMAX 10000

void EulerStep(double, double*, void(*) (double, double*, double*), double, int);
void RK2Step(double, double*, void(*) (double, double*, double*), double, int);
void RK4Step(double, double*, void(*) (double, double*, double*), double, int);
void dYdt(double, double*, double*);

double Residual(double);
double Analytic(double);

int main(){
  
  #if DEBUG == TRUE
  cout << "!! DEBUG TIME !!" << endl;
  #endif

  #if SESSION == 1
  int neq = 2;
  double Y[neq];
  double r = 0., s = 0.;
  int N = 1000; //numero di punti  
  double ri = 0., rf = 20.;
  double dr = (rf - ri)/(double)N;

  ofstream fdata;
  fdata.open("poisson.dat");
  for(int i = 0; i < 6; i++){
    Y[0] = 0.;
    Y[1] = s;
    for(int j = 0; j < N; j++){
      RK4Step(r,Y,dYdt,dr,neq);
      r += dr;    
      fdata << r << " " << Y[0] << endl;
    }
    fdata << endl << endl;
    r = 0.;
    s += .2;
  }
  fdata.close();

  s = 0.;
  double smax = 5.;
  double ds = (smax - s)/N;
  double res;

  fdata.open("residual.dat");
  do{
    res = Residual(s);
    fdata << s << " " << res << endl;
    s += ds;
  }while(s <= 5.);
  fdata.close();

  double root;
  double a = -10., b = 10.;
  double tol = 1.e-8;
  int k = 0;
  Bisection(Residual,a,b,tol,tol,k,root);
  cout << "Resiudal vanishes for s = " << root << endl;

  r = 0.;

  fdata.open("poisson_analytic.dat");
    Y[0] = 0.;
    Y[1] = root;
    for(int j = 0; j < N; j++){
      RK4Step(r,Y,dYdt,dr,neq);
      r += dr;    
      fdata << r << " " << Y[0] << " " << Y[1] << " " << Analytic(r) << endl;
    }
  fdata.close();
  #endif


  return 0;
}


void EulerStep(double t, double* Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){
  int k;
  double rhs[neq];
  
  RHSFunc(t,Y,rhs);
  for(k = 0; k < neq; k++)
    Y[k] += dt*rhs[k];
}

void RK2Step(double t, double* Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){
  double Y1[neq], k1[neq], k2[neq];

  RHSFunc(t,Y,k1);
  for(int i = 0; i < neq; i++){
    Y1[i] = Y[i] + 0.5*dt*k1[i];
  }

  RHSFunc(t+0.5*dt,Y1,k2);
  
  for(int i = 0; i < neq; i++){
    Y[i] += dt*k2[i];
  }
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
  double rho = exp(-t)/(8*M_PI);
  R[0] = v;
  R[1] = -4*M_PI*t*rho;
}

double Residual(double s){
  int neq = 2;
  double Y[neq];
  double r = 0.;
  int N = 1000; //numero di punti  
  double ri = 0., rf = 20.;
  double dr = (rf - ri)/(double)N;
  
  Y[0] = 0.;
  Y[1] = s;
  for(int j = 0; j < N; j++){
    RK4Step(r,Y,dYdt,dr,neq);
    r += dr;    
  }

  return Y[0] - 1.;
}

double Analytic(double r){
  return 1. - 0.5*(r + 2.)*exp(-r);
}
