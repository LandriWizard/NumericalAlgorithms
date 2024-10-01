#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define NSIZE 1000

void Average(double x[],int,double&,double&,double&);

int main(){
  int i;
  int n=NSIZE;
  double x[NSIZE];
  double xav, s2, s; //media, varianza, deviazione standard	
  
  srand48(time(NULL));
  
  for(i=0;i<n;i++){
    x[i] = drand48();
  }

  Average(x, n, xav, s2, s);
  cout << "<x> = " << xav << "; s^2 = " << s2 << "; s = " << s << endl;
  cout << "Expected : <x> = " << 0.5 << "; s^2 = " << 1.0/12.0 << "; s = " << 1.0/sqrt(12.0) << endl;
  return 0;
}

void Average(double x[], int n, double & xav, double & s2, double & s){
  double mu=0;
  for(int i=0;i<n;i++){
  mu += x[i];
  }
  xav = mu/n;
  s2=0;
  for(int i=0;i<n;i++){
    s2 += (x[i] - xav)*(x[i] - xav);
  }
  s2 /= n;
  s = sqrt(s2);
}
