#include <iostream>
using namespace std;

int Quotient(int,int,int,int);

int main(){
  int a, b, q, r;
  cout << "Dammi due interi (a,b)" << endl;
  cin >> a >> b;
  int flag = Quotient(a,b,q,r);
  if(flag!=0){
    cout << "!! Non si puo dividere per zero !!" << endl;
    return 1;
  }
  Quotient(a,b,q,r);
  cout << a << " / " << b << " ha parte reale " << q << " e resto " << r << endl;
  return 0;
}

int Quotient(int a, int b, int & q, int & r){

  if(b==0) return 1;
  q = a/b;
  r = a%b;
  return 0;

}
