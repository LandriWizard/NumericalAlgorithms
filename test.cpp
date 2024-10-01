#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#define TRUE 1
#define FALSE 0
#define DEBUG TRUE

#define ARRAYMAX 100

void EulerStep(double, double*, void(*) (double, double*, double*), double, int);
void dYdt(double, double*, double*);

int main(){

  cout << setiosflags(ios::scientific);
  cout << setprecision(4);

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
  fdata.open("test.dat");
  do{
    t = 0.;
    dt = (tmax - t)/(double)n;
    for(int i = 0; i < n; i++){
      EulerStep(t,Y,dYdt,dt,neq);
      t += dt;
      #if DEBUG == TRUE
      cout << "Y " << Y[0] << endl;
      #endif
    }
    errE = fabs(Y[0] - cos(t));
    #if DEBUG == TRUE
    cout << "t " << t << endl;
    cout << "errE " << errE << endl;
    #endif
    n *= 2;
  }while(n <= 16);

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


void dYdt(double t, double* Y, double* R){
  double x = Y[0];
  double y = Y[1];
  R[0] = y;
  R[1] = -x;
}

