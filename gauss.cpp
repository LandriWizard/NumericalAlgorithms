#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
using namespace std;

int main(){
  
  cout << setiosflags(ios::scientific);
  cout << setprecision(14);

  double sigma = 0.5; //sigma della gaussiana
  double xr; //ascissa casuale
  double yr; //ordinata casuale
  double f;

  srand48(time(NULL));  //seed del rand

//definisco il file che manderò a gnuplot
  ofstream fdata;
  fdata.open("gauss.dat");

  int n = 1.e5; //numero di coppie casuali che calcolo

  for(int k = 0; k < n; k++){
    xr = 10. * (drand48() - 0.5);
    yr = drand48();
    f = 1 / (sigma * sqrt(2*M_PI)) * exp(-0.5*xr*xr/(sigma*sigma));
    if(yr < f){
      fdata << xr << " " << yr << endl;
    }
  }
  fdata.close();

  return 0;
}
