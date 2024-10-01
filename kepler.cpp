#define _USE_MATH_DEFINES 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

#define FALSE 0
#define TRUE 1
#define DEBUG FALSE 

void RK4Step(double, double*, void(*) (double, double*, double*), double, int);
void dYdt(double, double*, double*);

int main(){
  
  #if DEBUG == TRUE
  cout << "!! DEBUG TIME !!" << endl;
  #endif

  cout << setiosflags(ios::scientific);
  cout << setprecision(4);

  int neq = 4, max_nsteps = 100000;
  int i;
  double dt = 0.01, t, r, v, T;
  double Y[neq];
  double vy_old; //variabile che mantiene in memoria un valore precedente
  int nMaxOrbits = 10;
  int nInversions = 0;
  double x0, y0;
  double vx0, vy0;
  double E0, E;

  ofstream fdata;
  fdata.open("kepler.dat");
//1. Assign initial conditions and calculate energy
  x0 = 4.;
  y0 = 0.;
  r = sqrt(x0*x0 + y0*y0);
  vx0 = 0.;
  vy0 = sqrt(0.3/r);
  v = sqrt(vx0*vx0 + vy0*vy0);
  E0 = 0.5*v*v - 1./r;

  Y[0] = x0, Y[1] = y0;  //posizioni iniziali (x,y)
  Y[2] = vx0, Y[3] = vy0; //velocità iniziali in x e y
/*while(1){  
    yP = Y[2];
    RK4Step(t,Y,dYdt,dt,neq);
    t += dt;    
    #if DEBUG == TRUE
    cout << nOrbits << "\t" << Y[0] << "\t" << Y[1] << "\t" << Y[2] << "\t" << Y[3] << endl;
    #endif
    fdata << t << " " << Y[0] << " " << Y[1] << endl; 
    if(yP*Y[2] <= 0)
      nInv++;
    if(nInv%2 == 0 &&  nInv != 0){
      nOrbits++;
      nInv = 0;
    }
    if(nOrbits == nMaxOrbits)
      break;
  }*/
//2. Integrate with RK4
  t = 0.;
  for(i = 1; i <= max_nsteps; i++){
    //dt = 0.5*r/v;
    dt = 0.05*r/v;
    vy_old = Y[3];
    RK4Step(t,Y,dYdt,dt,neq);
    t += dt;
    if(Y[3]*vy_old < 0.){
      nInversions++;
      cout << "RK Step = " << i << "; turning point #" << nInversions << endl;
      if(nInversions == 2*nMaxOrbits + 1) break;
    }
    r = sqrt(Y[0]*Y[0] + Y[1]*Y[1]);
    v = sqrt(Y[2]*Y[2] + Y[3]*Y[3]);
    #if DEBUG == TRUE
      cout << "r = " << r << ";\tv = " << v << endl;
    #endif
    E = 0.5*v*v - 1./r;
    fdata << t << " " << Y[0] << " " << Y[1] << " " << dt << " " << fabs(Y[0]*Y[3] - Y[1]*Y[2]) << " " << fabs(E/E0 - 1.) << endl;
  }
  fdata.close();
 
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
  double x = Y[0]; //posizione x
  double y = Y[1]; //posizione y
  double vx = Y[2]; //velocità x
  double vy = Y[3]; //velocità y
  double r2 = x*x + y*y; //raggio
  double r3 = r2*sqrt(r2);
  R[0] = vx;
  R[1] = vy;
  R[2] = -x/r3;
  R[3] = -y/r3;
}
