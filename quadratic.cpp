#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

int main(){
  
  double a=1, b, c, x1, x2; //coefficienti dell'equazione e soluzioni date dall'utente
  double x1c, x2c; //soluzioni calcolate per un controllo della precisione
  
  cout << "Dimmi le soluzioni dell'equazione quadratica che vuoi controllare" << endl;
  cin >> x1 >> x2;
  b = - ( x1 + x2 );
  c = x1 * x2;

  double delta = b * b - 4 * a * c;
  
  cout << "Calcolo con la formula risolutiva standard" << endl;

  x1c = 0.5 * ( - b - sqrt(delta) ) / a;//calcolo x1c
  x2c = 0.5 * ( - b + sqrt(delta) ) / a;//calcolo x2c
  
  //cout << "L'equazione e' " << a << "x^2 + " << b << "x + " << c << " = 0" << endl;
  cout << "x1 = " << x1 << "; x1c = " << x1c << endl;
  cout << "x2 = " << x2 << "; x2c = " << x2c << endl;

  if(x1c == x1 && x2c == x2)
    cout << "I risultati sono esatti" << endl;
  else
    cout << "I risultati sono errati" << endl;
  
  cout << "Calcolo razionalizzando i meno" << endl;

  double x1cc, x2cc;

  if(b >= 0){
    x1cc = 0.5 * ( - b - sqrt(delta) ) / a;
    x2cc = (2 * c) / ( - b - sqrt(delta) );
  }
  else{
    x1cc = (2 * c) / ( - b + sqrt(delta) );
    x2cc = 0.5 * ( - b + sqrt(delta) ) / a;
  }

  //cout << "L'equazione e' " << a << "x^2 + " << b << "x + " << c << " = 0" << endl;
  cout << "x1 = " << x1 << "; x1cc = " << x1cc << endl;
  cout << "x2 = " << x2 << "; x2cc = " << x2cc << endl;

  if(x1cc == x1 && x2cc == x2)
    cout << "I risultati sono esatti" << endl;
  else
    cout << "I risultati sono errati" << endl;

  return 0;
}
