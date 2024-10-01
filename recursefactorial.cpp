#include <iostream>
#include <climits>
#include <cfloat>
using namespace std;

int Factorial(int);
double Dfactorial(int);
//double Ifactorial(int, double);
//double Idfactorial(int, double);

int main(){
  int n; //numero di cui calcolo il fattoriale
  double f; //fattoriale intero
  double df; //fattoriale double
  cout << "Dimmi un intero di cui vuoi calcolare il fattoriale" << endl;
  cin >> n;
  if(n == 0 || n == 1){
    cout << n << "! = " << 1 << endl;
    return 1;
  }
  //f = (double) Factorial(n); //usato per il fattoriale inverso
  f = Factorial(n); //confrontato con int max
  df = Dfactorial(n);
  cout << "Usando gli interi" << endl;
  cout << n << "! = " << Factorial(n) << endl;
  cout << "Usando i double" << endl;
  cout << n << "! = " << Dfactorial(n) << endl;
  
  cout << "Check di validita' dei risultati" << endl;
  cout << "Per numeri interi" << endl;
  if(f < INT_MAX)
    cout << "Risultato giusto" << endl;
  else
    cout << "Numero troppo grande per calcolarne il fattoriale con gli int" << endl;
  //Ifactorial(n,f);

  cout << "Per double" << endl;
  if(df < DBL_MAX)
    cout << "Risultato giusto" << endl;
  else
    cout << "Numero troppo grande per calcolarne il fattoriale con i double" << endl;
  //Idfactorial(n,df);

  return 0;
}

int Factorial(int n){
  if(n > 1)
    return n*Factorial(n - 1);
  else
    return 1;
}

double Dfactorial(int n){
  if(n > 1)
    return n*Dfactorial(n - 1);
  else
    return 1;
}

/*double Ifactorial(int n, double f){
  if(f > 1){
    //cout << f << endl;
    f /= n;
    return Ifactorial(n-1,f);
  }
  else if(f == 1){
    cout << "Risultato giusto" << endl;
    return 1;
  }
  else{
    cout << "Numero troppo grande per calcolarne il fattoriale con gli interi" << endl;
    return 1;
  }
}

double Idfactorial(int n, double df){
  if(df > 1){
    //cout << df << endl;
    df /= n;
    return Idfactorial(n-1,df);
  }
  else if(df == 1){
    cout << "Risultato giusto" << endl;
    return 1;
  }
  else{
    cout << "Numero troppo grande per calcolarne il fattoriale con i double" << endl;
    return 1;
  }
}*/
