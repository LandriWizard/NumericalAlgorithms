#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
using namespace std;

int main(){
  
  cout << setiosflags(ios::scientific);
  cout << setprecision(14);

  long int N = 1000;

  //double r[N];
  
  double r, k1, k2;

  srand48(time(NULL));  

  ofstream fdata;
  fdata.open("uniformity.dat");
  for(int i = 0; i < N; i++){
    r = drand48();
    fdata << r << " " << i << endl;
  }
  fdata.close();

  fdata.open("prn_uniformity.dat");
  for(int N = 4; N <= 1.e7; N *= 2){
    k1 = 0;
    k2 = 0;
    for(int i = 0; i < N; i++){
      r = drand48();
      k1 += r;    //momento 1
      k2 += r*r;  //momento 2
    }
  k1 /= (double)N;
  k2 /= (double)N;
  fdata << N << " " << fabs(k1 - 0.5) << " " << fabs(k2 - 1./3.) <<endl;
  cout << N << " " << fabs(k1 - 0.5) << " " << fabs(k2 - 1./3.) <<endl;
  }
  fdata.close();

  return 0;
}
