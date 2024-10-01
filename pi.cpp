#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
#include <climits>
using namespace std;

int main(){
  
  cout << setiosflags(ios::scientific);
  cout << setprecision(14);

  double tol = 1.e-4; //tolleranza
  double x; //ascissa casuale
  double y; //ordinata casuale
  double sigma;
  int cnt = 0; //conta i punti dentro al cerchio
  double pi; //integrale che voglio approssimare
  double pi2;
  double err = 1.; //errore dell'approssimazione
  int n;
  
  double iseed = time(NULL);
  srand48(iseed);  //seed del rand

//definisco il file che manderò a gnuplot
  ofstream fdata;

//simple integration  using fixed number of points

  iseed = time(NULL);
  srand48(iseed);

  n = 600;
  cnt = 0;
  for(int i = 0; i < n; i++){
    x = (drand48()-0.5)*2;
    y = (drand48()-0.5)*2;
    if ( (x*x * y*y) <= 1. ) cnt++;
    //cout << x << "\t\t" << y << endl;
  }
  pi = 4.*(double)cnt/n;
  cout << "\n1. Fix value of N and compute pi" << endl;
  cout << "     n = " << n
       << ";   pi = " << pi
       << ";  err = " << fabs(pi/(M_PI-1.)) << endl;

//2. Now determine how many tiral steps ae necessary to reach specified tolerance tol
  fdata.open("pi.dat");

  cout << "Now computing how many steps are necessary to reach the tolerance of " << tol << endl;
  cnt = n = 0;

  do{
    n++;
    x = (drand48()-0.5)*2;
    y = (drand48()-0.5)*2;
    if ( (x*x + y*y) <= 1. ) cnt++;
    pi = 4.*(double)cnt/n;
    err = fabs(pi/(M_PI-1.));
    fdata << x << " " << y << endl; 
    cout << n << "\t\t " << pi << "\t\t" << err << endl;
    //if(n == 20) break;
  }while(err > tol);

  fdata.close();

//3. Now determine the error |I-pi| as a function of N by doubling N

  double Asq = 4.;
  fdata.open("pi_err.dat");
  cout << "\n3. Now computing error" << endl;
  n = 4;
  for(int i = 1; i < 24; i++){
    
    if(n > INT_MAX) {
      cout << "! Integer overflow" << endl;
      return 1;
    }
    srand48(time(NULL));
    cout << "n = " << n << " (i = " << i << ")" << endl;
    cnt = 0;

    pi = pi2 = 0.;
    for(int k = 0; k < n; k++){
      x = (drand48() - 0.5) * 2;
      y = (drand48() - 0.5) * 2;
      if((x*x + y*y) <= 1){
        cnt++;
        pi += Asq;
        pi2 += Asq*Asq;
      }
    }
    pi /= n;
    pi2 /= n;

    double varDist = n/(n-1.) * (pi2 - pi*pi);
    sigma = sqrt(varDist/n);
    
    err = fabs(pi/M_PI-1.);
    fdata << n << " " << err << " " << sigma << endl;
    n *= 2;
  }
  fdata.close();

  return 0;
}
