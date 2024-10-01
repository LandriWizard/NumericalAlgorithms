#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

int main(){

  cout << setiosflags(ios::scientific);
  cout << setprecision(8);

  /*double x0, x1, x2;
  double h = 0.5;
  double FD, BD, CD;
  double errFD, errBD, errCD;
  double flag;
  double tol = 1.e-7;

  ofstream fdata;
  fdata.open("derivative.dat");
  do{
    x0 = sin(1.);
    x1 = sin(1.- h);
    x2 = sin(1. + h);
    FD = (x2 - x1)/h;
    BD = (x1 - x0)/h;
    CD = 0.5*(x2 - x0)/h;
    errFD = fabs(FD - cos(1.));
    errBD = fabs(BD - cos(1.));
    errCD = fabs(CD - cos(1.));
    flag = 1./h;
    fdata << flag << " " << errFD << " " << errBD << " " << errCD << " " << endl;
    h /= 2.;
  }while(errCD > tol);
  fdata.close(); 
  return 0;*/

  double x = 1.;
  double truedev = cos(x);
  double fd, bd, cd;
  double tol = 1.e-7;
  double h = 0.5, flag;
  ofstream fdata;
  fdata.open("derivative.dat");
  do{
    fd = (sin(x + h) - sin(x))/h;
    bd = (sin(x) - sin(x - h))/h;
    cd = 0.5*(sin(x + h) - sin(x - h))/h;
    flag = 1./h;
    fdata << flag << " " << fabs(fd - truedev) << " " << fabs(bd - truedev) << " " << fabs(cd - truedev) << endl;
    h *= 0.5;
  }while(fabs(cd-truedev) > tol);
  fdata.close();
  return 0;
}
