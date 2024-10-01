#include <iostream>
using namespace std;

int main(){

  double dnum = 1.;
  double deps = 1.;

  while(dnum + deps != dnum){
    deps /= 10;
  }
  
  cout << "Precisione di macchina double = " << deps << endl;

  float fnum = 1.;
  float feps = 1.;

  while(fnum + feps != fnum){
    feps /= 10;
  }
  
  cout << "Precisione di macchina float = " << feps << endl;
  
  return 0;
}
