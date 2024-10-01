#include "my_header.h"

double PolyLeg(double x, int n){
  if(n == 0)
    return 1.;
  if(n == 1)
    return x;
  else
    return ((2*n-1)*x*PolyLeg(x,n-1)-(n-1)*PolyLeg(x,n-2))/n;
}

double DPolyLeg(double x, int n){
  if(n == 0)
    return 1.;
  else
    return (n/(x*x - 1)*(x*PolyLeg(x,n) - PolyLeg(x,n-1)));
}
