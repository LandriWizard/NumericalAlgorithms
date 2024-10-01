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

static double g_k = 1.;

void RK4Step(double, double*, void(*) (double, double*, double*), double, int);
void dYdt(double, double*, double*);

double Residual(double);

int main(){
  
  #if DEBUG == TRUE
  cout << "!! DEBUG TIME !!" << endl;
  #endif

  #if SESSION == 1
  double a = 1.; 
  double b = 20.;
  int nInters = 41;
  double xL[nInters], xR[nInters];
  double root;
  int iterations;
  int nRoots;
  Bracket(Residual,a,b,xL,xR,nInters,nRoots);
  int j = 0;
  double rootsArr[nRoots];
  double code;
  for(int i = 0; i < nRoots; i++){
    code = Bisection(Residual,xL[i],xR[i],1.e-8,1.e-8,iterations,root);
    if(code == 0){
      rootsArr[j] = root;
      j++;
    }
  }
  double start = 0., end = 1.;
  int N = 100;
  double h = fabs(end - start)/(double)N;
  int neq = 2;
  double Y[neq];
  double t;
  ofstream fdata;
  fdata.open("wave.dat");
  for(int i = 0; i < nRoots; i++){
    t = 0;
    Y[0] = 0., Y[1] = 1.;
    g_k = rootsArr[i];
    for(int j = 0; j < N; j++){
      fdata << t << " " << Y[0] << " " << Y[1] << endl;
      RK4Step(t,Y,dYdt,h,neq);
      t += h;
    }
    fdata << t << " " << Y[0] << " " << Y[1] << endl;
    fdata << endl << endl;
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
  R[1] = -g_k*g_k*x;
}

double Residual(double k){
  g_k = k;
  double start = 0.;
  double end = 1.;
  double t = 0.;
  int N = 100; //numero di punti  
  double h = fabs(end - start)/(double)N;
  int neq = 2;
  double Y[neq];
  Y[0] = 0.;
  Y[1] = 1.;
  for(int j = 0; j < N; j++){
    RK4Step(t,Y,dYdt,h,neq);
    t += h;    
  }

  return Y[0];
}
