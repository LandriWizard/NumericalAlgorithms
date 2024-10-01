#include <iostream>
using namespace std;

int main(){
  
  //Esercizio con i float
  
  float a,b; //inizializzo i due numeri
  cout << "Dammi due numeri reali" << endl;
  cout << "a = ";
  cin >> a; //a da tastiera
  cout << "\nb = ";
  cin >> b; //b da tastiera

  float sumf = a + b;
  cout << a << " + " << b << " = " << sumf << endl;
  float diff = a - b;
  cout << a << " - " << b << " = " << diff << endl;
  float mulf = a * b;
  cout << a << " * " << b << " = " << mulf << endl;
  float divf = a / b;
  cout << a << " / " << b << " = " << divf << endl;
  
  // esercizio con gli int

  int c,d; //inizializzo i numeri
  cout << "Dammi due numeri interi" << endl;
  cout << "c = ";
  cin >> c; //a da tastiera
  cout << "\nd = ";
  cin >> d; //b da tastiera

  int sumi = c + d;
  cout << c << " + " << d << " = " << sumi << endl;
  int difi = c - d;
  cout << c << " - " << d << " = " << difi << endl;
  int muli = c * d;
  cout << c << " * " << d << " = " << muli << endl;
  float divi = (float)c / (float)d; 
  cout << c << " / " << d << " = " << divi << endl;
}
