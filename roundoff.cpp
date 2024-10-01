#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

int main(){
  
  cout << "Example #1: compute sqrt(x^2 + 1) - x, for large values of x" << endl;

  for(float i = 1.e4; i <= 1.e10 ; i *= 10){
    float fx11 = sqrt( i * i + 1.) - i; //funzione esatta
    float fx21 = 1. / (sqrt(i*i + 1.) + i); //funzione razionalizzata
    float fxt1 = 0.5 * (1. / i); //serie di Taylor ordine 2
    cout << "x = " << i << "; fx1 = " << fx11 << "; fx2 = " << fx21 << "; f(taylor) = " << fxt1 << endl;  
  }

  cout << "Example #2: compute 1 - cos(x), for small values of x" << endl;
  for(float j = 1.e-1; j > 1.e-9; j /= 10){
    float fx12 = 1. - cos(j); //funzione esatta
    float fx22 = ( sin(j) * sin(j) ) / ( 1. + cos(j)); //funzione razionalizzata
    float fxt2 = 0.5 * j * j; //serie di Taylor ordine 2
    cout << "x = " << j << "; fx1 = " << fx12 << "; fx2 = " << fx22 << "; f(taylor) = " << fxt2 << endl;  
  }

  return 0;

}
