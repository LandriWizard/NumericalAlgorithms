#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
#include <climits>
#include "my_header.h"
using namespace std;

#define DIMMAX 24
#define ARRAYMAX 256

#define SESSION 5 

int g_LengendreN = 5;

double Func(double);
double DFunc(double);
double Bisection(double (*F) (double), double, double, double);
double FalsePos(double (*F) (double), double, double, double);
double Secant(double (*F) (double), double, double, double);
double Newton(double (*F) (double), double (*G) (double), double, double, double, int&, double&);
void Bracket(double (*F) (double), double, double, double*, double*, int, int&);

int main(){
  
  cout << setiosflags(ios::scientific);
  cout << setprecision(4);
  
  double a, b, tol;

  #if SESSION == 1 
  cout << "\n--------------- SESSION 1 --------------\n" << endl;

  a = -1., b = 1.; //limiti dell'intervallo in cui cerco lo zero
  tol = 1.e-7; //tolleranza sullo zero
  
  cout << "Finding zeros with Bisection algorithm\n" << endl;
  Bisection(Func,a,b,tol);

  cout << "Using False position algorithm\n" << endl;
  FalsePos(Func,a,b,tol);

  cout << "Using Secant algorithm\n" << endl;
  Secant(Func,a,b,tol);

  cout << "Using Newton algorithm\n" << endl;   
  Newton(Func,DFunc,a,b,tol);
  #endif  

  #if SESSION == 3
  cout << "\n--------------- SESSION 3 --------------\n" << endl;

  a = -5., b = 0.; //nuovo intervallo fer la funzione polinomiale
  tol = 1.e-8;     //nuova tolleranza per la polinomiale

  cout << "Finding zeros with Bisection algorithm\n" << endl;
  Bisection(Func,a,b,tol);

  cout << "Using False position algorithm\n" << endl;
  FalsePos(Func,a,b,tol);

  cout << "Using Secant algorithm\n" << endl;
  Secant(Func,a,b,tol);

  cout << "Using Newton algorithm\n" << endl;   
  Newton(Func,DFunc,a,b,tol);
  #endif

  #if SESSION == 4
  cout << "\n--------------- SESSION 4 --------------\n" << endl;

  a = 0.; b = 2.;
  tol = 1.e-7;

  cout << "Finding with Bisection algorithm\n" << endl;
  Bisection(Func,a,b,tol);
  cout << "Using False position algorithm\n" << endl;
  FalsePos(Func,a,b,tol);

  cout << "Using Secant algorithm\n" << endl;
  Secant(Func,a,b,tol);

  cout << "Using Newton algorithm\n" << endl;   
  Newton(Func,DFunc,a,b,tol);
  #endif

  #if SESSION == 5
  cout << "\n--------------- SESSION 5 -------------\n" << endl;

  a = -10.0, b = 10.0;
  double xtol = 1.e-8, root = 0.0;
  int n = 10, nroots, r, ntry;
  double xL[ARRAYMAX], xR[ARRAYMAX];

  Bracket(Func, a, b, xL, xR, n, nroots);

  for(int i = 0; i < nroots; i++){
    cout << "#" << i+1 << "; xL = " << xL[i] << "; xR = " << xR[i]<< endl;
  }
 
  for(int i = 0; i < nroots; i++){
    int r = Newton(Func, DFunc, xL[i], xR[i], xtol, ntry, root);
    cout << root << endl;
  }
  #endif

  #if SESSION == 6

  #endif

  return 0;
  
}


#if SESSION == 1
double Func(double x){
  return exp(-x) - x;
}

double DFunc(double x){
  return -exp(-x) - 1.;
}
#endif

#if SESSION == 3
double Func(double x){
  return 5. + x*(1. + x*(-3. + x));
}

double DFunc(double x){
  return 1. + x*(-6. + 3*x); 
}
#endif 

#if SESSION == 4
double Func(double x){
  return exp(1/(x + 0.5)) - (3. + 2.*x)/(1 + x);
}

double DFunc(double x){
  return -exp(1./(x + 0.5))/(x + 0.5) + 1./(1. + 2.*x + x*x);
}
#endif

#if SESSION == 5
double Func(double x){
  return sin(x) - (1./3. + x/5. + x*x/100.);
}

double DFunc(double x){
  return cos(x) -  (x/50. + 1./5.);
}
#endif

#if SESSION == 6
  double Func(double x){
    double px;
    for(int i = 0; i < g_LegendreN; i++){
      if(i = 0)
        px = 1.;
      else if(i = 1)
        px = x;
      else{
        

      }        
    }
  }
#endif
