#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main(){

  double s; //numero di cui calcoliamo la radice quadrata
  double x; //numero che stimiamo sia la radice quadrata di s
  double err; //numero che ci dice la differenza tra due interazioni successive
  double xp; //valore in cui memorizzo l'iterazione precedente
  int i; //counter delle iterazioni

  cout << "Dammi un numero di cui calcolare la radice quadrata con il metodo di Erone:" << endl;
  cin >> s;
  cout << "Dammi una stima della radice quadrata:" << endl;
  cin >> x;

  for(i = 0; 1; i++){
    xp = x; //fisso il valore di x nell'iterazione precedente
    x = 0.5 * (x + s/x); //aggiorno x a uno piu' vicino al risultato esatto
    err = fabs(x - xp); //calcolo l'errore dell'iterazione precedente rispetto a quella attuale
    cout << scientific << setw(22) << setprecision(14) << "Iterazione #" << i+1 << "; x = " << x << "; err =" << err << endl;
    if( err < 1.e-15 ) break; 
  }

  cout << scientific << setw(22) << setprecision(14) << "La radice quadrata di " << s << " e' :" << x << endl;
  cout << scientific << setw(22) << setprecision(14) << "Il valore esatto e': " << sqrt(s) << endl;

  return 0; 
}
