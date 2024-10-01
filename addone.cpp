#include <iostream>
using namespace std;

float addone(float);

int main(){
  float a;
  cout << "Dammi un numero reale" << endl;
  cin >> a;
  cout << a << " + 1 = " << addone(a) << endl;
  return 0;
}

float addone(float x){
  return x+=1;
}
