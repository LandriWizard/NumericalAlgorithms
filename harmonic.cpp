#define _USE_MATH_DEFINES 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

#define FALSE 0
#define TRUE 1
#define DEBUG TRUE

#define ARRAYMAX 10000

void EulerStep(double, double*, void(*) (double, double*, double*), double, int);
void RK2Step(double, double*, void(*) (double, double*, double*), double, int);
void RK4Step(double, double*, void(*) (double, double*, double*), double, int);
void PosVerlet(double, double*, double*, void(*) (double, double*, double*, int), double, int);

void Acceleration(double, double*, double*, int);
void dYdt(double, double*, double*);

int main(){
  
  #if DEBUG == TRUE
  cout << "!! DEBUG TIME !!" << endl;
  #endif

  cout << setiosflags(ios::scientific);
  cout << setprecision(4);

  int N = 1;
  double T, t = 0.;
  T = 2*M_PI;
  double dt = 0.02*T;

//POSIITION VERLET
  double X[N]; //posizioni
  double V[N]; //velocità
  X[0] = 1.;
  V[0] = 0.;
  int nMaxOscil = 10;
  double E;
  ofstream fdata;
  fdata.open("harmonic.dat");
  double tenthT = 100*T;
  for(int i; i < tenthT; i++){
    PosVerlet(t,X,V,Acceleration,dt,N);
    t += dt;
    E = .5*(V[0]*V[0] + X[0]*X[0]);
    fdata << t << " " << X[0] << " " << V[0] << " " << E << endl;
  }
  fdata << endl << endl;
  
//RK2
  int neq = 2;
  t = 0.;
  double Y[neq];
  Y[0] = 1., Y[1] = 0.;
  for(int i = 0; i < tenthT; i++){
    RK2Step(t,Y,dYdt,dt,neq);
    t += dt;
    E = .5*(Y[0]*Y[0] + Y[1]*Y[1]);
    fdata << t << " " << Y[0] << " " << Y[1] << " " << E << endl;
  }
  fdata.close();

 
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
  double y0 = Y[0];
  double y1 = Y[1];
  R[0] = y1;
  R[1] = -y0;
}

void Acceleration(double t, double* X, double* A, int N){
  double x0 = X[0];
  A[0] = -x0;
}

void PosVerlet(double t, double* X, double* V, void (*RHSFunc)(double, double*, double*, int), double dt, int N){
  double A[N];
  for(int i = 0; i < N; i++)
    X[i] += .5*dt*V[i];
  RHSFunc(t,X,A,N);
  for(int i = 0; i < N; i++)
    V[i] += dt*A[i];
  for(int i = 0; i < N; i++)
    X[i] += .5*dt*V[i];
}
