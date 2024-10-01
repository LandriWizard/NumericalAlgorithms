#define _USE_MATH_DEFINES 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

#define FALSE 0
#define TRUE 1
#define DEBUG TRUE
#define SESSION 3

#define ARRAYMAX 10000

void EulerStep(double, double*, void(*) (double, double*, double*), double, int);
void RK2Step(double, double*, void(*) (double, double*, double*), double, int);
void RK4Step(double, double*, void(*) (double, double*, double*), double, int);
void dYdt(double, double*, double*);

int main(){
  
  #if DEBUG == TRUE
  cout << "!! DEBUG TIME !!" << endl;
  #endif

  cout << setiosflags(ios::scientific);
  cout << setprecision(4);

  #if SESSION == 1

  int neq = 1;
  double Y[ARRAYMAX];
  Y[0] = 1.;
  double t = 0., tmax = 3.;
  double dt = 0.2;
  double n;

  n = round((tmax - t)/dt);
  #if DEBUG == TRUE
  cout << n << endl;
  #endif

  #if DEBUG == TRUE
  cout << endl;
  cout << "x\t" << "\ty(x)\t" << "\tabs_err\t" << "\trel_err\t" << endl;
  cout << t << "\t" << Y[0] << "\t" << fabs(Y[0] - exp(-0.5*t*t)) << "\t" << fabs((Y[0] - exp(-0.5*t*t))) << "\t" << endl;
  #endif

  ofstream fdata;
  fdata << setiosflags(ios::scientific);
  fdata << setprecision(4);
  fdata.open("ode1_Euler.dat");
  fdata << t << " " << Y[0] << " " << fabs(Y[0] - exp(-0.5*t*t)) << " " << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << endl;
  for(int i = 0; i < n; i++){
    EulerStep(t,Y,dYdt,dt,neq);
    t += dt;
    #if DEBUG == TRUE
    cout << t << "\t" << Y[0] << "\t" << fabs(Y[0] - exp(-0.5*t*t)) << "\t" << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << "\t" << endl;
    #endif
    fdata << t << " " << Y[0] << " " << fabs(Y[0] - exp(-0.5*t*t)) << " " << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << endl;
  }
  fdata.close();  
  #if DEBUG == TRUE 
  cout << endl;
  #endif

//-------------------------------------------------------------------RK2---------------------------------------------------//
  
  #if DEBUG == TRUE
  cout << "!! RUNGE KUTTA 2ND ORDER !!" << endl;
  #endif

  Y[0] = 1.;
  t = 0.;

  #if DEBUG == TRUE 
  cout << endl;
  cout << "x\t" << "\ty(x)\t" << "\tabs_err\t" << "\trel_err\t" << endl;
  cout << t << "\t" << Y[0] << "\t" << fabs(Y[0] - exp(-0.5*t*t)) << "\t" << fabs((Y[0] - exp(-0.5*t*t))) << "\t" << endl;
  #endif

  fdata.open("ode1_RK2.dat");
  fdata << t << " " << Y[0] << " " << fabs(Y[0] - exp(-0.5*t*t)) << " " << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << endl;
  for(int i = 0; i < n; i++){
    RK2Step(t,Y,dYdt,dt,neq);
    t += dt;    
    #if DEBUG == TRUE
    cout << t << "\t" << Y[0] << "\t" << fabs(Y[0] - exp(-0.5*t*t)) << "\t" << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << "\t" << endl;
    #endif
    fdata << t << " " << Y[0] << " " << fabs(Y[0] - exp(-0.5*t*t)) << " " << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << endl;
  }
  fdata.close();
 
//-----------------------------------------------------------------------RK4-----------------------------------------------------------//

  #if DEBUG == TRUE
  cout << "!! RUNGE KUTTA 4ND ORDER !!" << endl;
  #endif

  Y[0] = 1.;
  t = 0.;

  #if DEBUG == TRUE 
  cout << endl;
  cout << "x\t" << "\ty(x)\t" << "\tabs_err\t" << "\trel_err\t" << endl;
  cout << t << "\t" << Y[0] << "\t" << fabs(Y[0] - exp(-0.5*t*t)) << "\t" << fabs((Y[0] - exp(-0.5*t*t))) << "\t" << endl;
  #endif

  fdata.open("ode1_RK4.dat");
  fdata << t << " " << Y[0] << " " << fabs(Y[0] - exp(-0.5*t*t)) << " " << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << endl;
  for(int i = 0; i < n; i++){
    RK4Step(t,Y,dYdt,dt,neq);
    t += dt;    
    #if DEBUG == TRUE
    cout << t << "\t" << Y[0] << "\t" << fabs(Y[0] - exp(-0.5*t*t)) << "\t" << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << "\t" << endl;
    #endif
    fdata << t << " " << Y[0] << " " << fabs(Y[0] - exp(-0.5*t*t)) << " " << fabs(Y[0] - exp(-0.5*t*t))/exp(-0.5*t*t) << endl;
  }
  fdata.close();
  #endif

  #if SESSION == 2
  int neq = 2;
  double Y[ARRAYMAX];
  Y[0] = 1., Y[1] = 0.;
  double t = 0., tmax = 20.*M_PI;
  int n = 200;
  double dt;
  dt = (tmax - t)/(double)n;
  
  #if DEBUG == TRUE
  cout << dt << endl;
  #endif
 
  ofstream fdata;
  fdata << setiosflags(ios::scientific);
  fdata << setprecision(4);
  fdata.open("ode2_Euler.dat");
  fdata << t << " " << Y[0] << " " << Y[1] << endl;
  for(int i = 0; i < n; i++){
    EulerStep(t,Y,dYdt,dt,neq);
    t += dt;
    #if DEBUG == TRUE
    cout << Y[0] << "\t" << Y[1] << endl;
    #endif
    fdata << t << " " << Y[0] << " " << Y[1] <<  endl;
  }
  fdata.close();  
  #if DEBUG == TRUE 
  cout << endl;
  #endif
//-------------------------------------------------------------------RK2---------------------------------------------------//
  
  #if DEBUG == TRUE
  cout << "!! RUNGE KUTTA 2ND ORDER !!" << endl;
  #endif

  Y[0] = 1., Y[1] = 0.;
  t = 0.;

  #if DEBUG == TRUE 
  cout << endl;
  cout << "x\t" << "\ty(x)\t" << endl;
  cout << Y[0] << "\t" << Y[1] << endl;
  #endif

  fdata.open("ode2_RK2.dat");
  fdata << t << " " << Y[0] << " " << Y[1] << endl;
  for(int i = 0; i < n; i++){
    RK2Step(t,Y,dYdt,dt,neq);
    t += dt;    
    #if DEBUG == TRUE
    cout << Y[0] << "\t" << Y[1] << endl;
    #endif
    fdata << t << " " << Y[0] << " " << Y[1] << endl;
  }
  fdata.close();
 //-----------------------------------------------------------------------RK4-----------------------------------------------------------//

  #if DEBUG == TRUE
  cout << "!! RUNGE KUTTA 4ND ORDER !!" << endl;
  #endif

  Y[0] = 1., Y[1] = 0.;
  t = 0.;

  #if DEBUG == TRUE 
  cout << endl;
  cout << "x\t" << "\ty(x)\t" << endl;
  cout << t << "\t" << Y[0] << endl;
  #endif

  fdata.open("ode2_RK4.dat");
  fdata << t << " " << Y[0] << " " << Y[1] << endl;
  for(int i = 0; i < n; i++){
    RK4Step(t,Y,dYdt,dt,neq);
    t += dt;    
    #if DEBUG == TRUE
    cout << Y[0] << "\t" << Y[1] << endl;
    #endif
    fdata << t << " " << Y[0] << " " << Y[1] << endl;
  }
  fdata.close();
  #endif
 
  #if SESSION == 3
  int neq = 2;
  double Y[ARRAYMAX];
  Y[0] = 1., Y[1] = 0.;
  double t = 0., tmax = 3.;
  int n = 4;
  double dt;
  double errE, errRK2, errRK4;
  ofstream fdata;
  fdata << setiosflags(ios::scientific);
  fdata << setprecision(4);
  fdata.open("ode2_Conv.dat");
  do{
    t = 0.;
    dt = (tmax - t)/(double)n;
    #if DEBUG == TRUE
    cout << "n " << n << endl;
    cout << "dt " << dt << endl;
    #endif
    Y[0] = 1., Y[1] = 0.;
    for(int i = 0; i < n; i++){
      EulerStep(t,Y,dYdt,dt,neq);
      t += dt;
    }
    errE = fabs(Y[0] - cos(t));
    #if DEBUG == TRUE
    cout << "EULER" << endl;
    cout << "dt " << dt << endl;
    cout << "errE " << errE << endl;
    #endif
    Y[0] = 1., Y[1] = 0.;
    t = 0.;
    for(int i = 0; i < n; i++){
      RK2Step(t,Y,dYdt,dt,neq);
      t += dt;    
    }
    errRK2 = fabs(Y[0] - cos(t));
    #if DEBUG == TRUE
    cout << "RK2" << endl;
    cout << "dt " << dt << endl;
    cout << "errRK2 " << errRK2 << endl;
    #endif
    Y[0] = 1., Y[1] = 0.;
    t = 0.;
    for(int i = 0; i < n; i++){
      RK4Step(t,Y,dYdt,dt,neq);
      t += dt;    
    }
    errRK4 = fabs(Y[0] - cos(t));
    #if DEBUG == TRUE
    cout << "RK4" << endl;
    cout << "dt " << dt << endl;
    cout << "errRK4 " << errRK4 << endl;
    cout << endl;
    #endif
    fdata << dt << " " << errE << " " << errRK2 << " " << errRK4 << endl;

    n *= 2;
  }while(n <= 2048);
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

#if SESSION == 1
void dYdt(double t, double* Y, double* R){
  double y = Y[0];
  R[0] = -t*y;
}
#endif

#if SESSION == 2
void dYdt(double t, double* Y, double* R){
  double x = Y[0];
  double y = Y[1];
  R[0] = y;
  R[1] = -x;
}
#endif

#if SESSION == 3
void dYdt(double t, double* Y, double* R){
  double x = Y[0];
  double y = Y[1];
  R[0] = y;
  R[1] = -x;
}
#endif

