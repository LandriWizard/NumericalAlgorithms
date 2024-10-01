#include <iostream>
using namespace std;

int sum(int, int);

int main(){

  int a,b;
  cout << "Dammi due interi (a,b)" << endl;
  cin >> a >> b;
  cout << a << " + " << b << " = " << sum(a,b) << endl;

  return 0;
}

int sum(int x, int y){
  return x+y;
}
