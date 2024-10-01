#include "my_header.h"

#define NGMAX 4

double GaussRule(double (*F)(double), double a, double b, int N, int Ng){
  
  double sum = 0.;
  double sumk;
  double x[NGMAX]; //array con le ascisse
  double w[NGMAX]; //array con i pesi
  double xv1;
  double xv2;
  double wv1;
  double wv2;
  double h = (b - a)/N;
  double x0;
  double x1;
  double deltax;
  double avx;
  
//riempio gli array in base al numero di punti di Gauss
  if(Ng == 2){
    xv1 = 1./sqrt(3.);
    x[0] = -xv1; x[1] = xv1;
    w[0] = 1; w[1] = 1;
  }
  else if(Ng == 3){
    xv1 = sqrt(3./5.);
    x[0] = -xv1; x[1] = 0; x[2] = xv1;
    w[0] = 5./9.; w[1] = 8./9.; w[2] = 5./9.;
  }
  else if(Ng == 4){
    xv1 = 0.339981043584856;
    xv2 = 0.861136311594053;
    x[0] = -xv1; x[1] = xv1; x[2] = -xv2; x[3] = xv2;
    wv1 = 0.652145154862546;
    wv2 = 0.347854845137454;
    w[0] = wv1; w[1] = wv1; w[2] = wv2; w[3] = wv2;
  }

  for(int i = 0; i < N; i++){
    x0 = a + i*h;
    x1 = a + (i+1)*h;
    deltax = 0.5*(x1 - x0);
    avx = 0.5*(x1 + x0);
    sumk = 0.;
    for(int k = 0; k < Ng; k++){
      sumk += w[k]*F(deltax*x[k] + avx);
    }
    sumk *= deltax;
    sum += sumk;
  }

  return sum;
}
