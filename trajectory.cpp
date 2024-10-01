#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

double Func(double);

#define FUNCTION 1 

int main(){

  cout << setiosflags(ios::scientific);
  cout << setprecision(5);

  double alpha = 10;
  double dt = 0.01*alpha;
  double v;
  double a;
  double t = 0;
  double fp;
  double f0;
  double fm;

  ofstream fdata;
  fdata.open("trajectory.dat");
  do{
    fp = Func(t + dt);
    f0 = Func(t);
    fm = Func(t - dt);    
    v = (fp - fm)/(2.*dt);
    a = (fp - 2*f0 + fm)/(dt*dt);
    fdata << t << " " << f0 << " " << v << " " << a << endl;
    t += dt;
  }while(t <= 10); 

  return 0;
}

#if FUNCTION == 1
double Func(double t){
  double a = 10;
  return a*t*t - t*t*t*(1. - exp((-a*a)/fabs(t)));
}
#endif

#if FUNCTION == 2
double Func(double t){
  double a = 10;
  return a*t*t - t*t*t*(1. - exp(fabs((-a*a)/t)));
}
#endif
