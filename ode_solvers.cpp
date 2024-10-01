#include "my_header.h"

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

