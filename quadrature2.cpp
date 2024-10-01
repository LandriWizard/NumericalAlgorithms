#include <iostream>
#include <iomanip>
#include <cmath>
#include "my_header.h"
using namespace std;

#define NGMAX 32

double Func(double);
double Root(double);
double Poly(double);
double RectangularRule(double (*F) (double), double, double, int);
double TrapezioidalRule(double (*F) (double), double, double, int);
double SimpsonRule(double (*F) (double), double, double, int);

int main(){
  double a = 0., b = 1.; //limiti dell'intervallo
  int N = 2; //numero di intervalli in cui divido l'intervallo originale
  double exact = 6.32120558813e-01;
  cout << "Exact:        " << scientific << setw(22) << setprecision(14) << exact << ";" << endl;
  double tol = 1.e-5;
  double sum = 0.;
  double sum0 = 0.;

//Calclolo con l'approssimazione rettangolare
  double err = 1.;
 
  while(err > tol){
    sum0 = sum;
    sum = RectangularRule(Func, a, b, N);
    err = fabs(sum - sum0);
    N *= 2;
  }
    cout << "Rectangular:  " << scientific << setw(22) << setprecision(14) << sum << "; N = " << N/2 << endl;

//Calclolo con l'approssimazione rettangolare
  N = 2;
  err = 1.;
  
  while(err > tol){
    sum0= sum;
    sum = TrapezioidalRule(Func, a, b, N); 
    err = fabs(sum - sum0);
    N *= 2;
  }
  cout << "Trapezioidal: " << scientific << setw(22) << setprecision(14) << sum << "; N = " << N/2 << endl; 

//Calcolo con l'approssimazione di Simpson
  N = 2;
  err = 1.;

  while(err > tol){
    sum0 = sum;
    sum = SimpsonRule(Func, a, b, N);
    err = fabs(sum - sum0);
    N *= 2;
  }
  cout << "Simpson:      " << scientific << setw(22) << setprecision(14) << sum << "; N = " << N/2 << endl;

//Adesso valuto la funzione Root con solo Simpson e Gauss

  cout << "\nCalcolo una funzione diversa" << endl;

//Calcolo con l'approssimazione di Simpson
  double c=0.;
  double d=3.;  
  N = 2;
  err = 1.;
  
  while(err > tol){
    sum0 = sum;
    sum = SimpsonRule(Root, c, d, N); 
    err = fabs(sum - sum0);
    N *= 2;
  }
  cout << "Simpson:      " << scientific << setw(22) << setprecision(14) << sum << ";" << endl;

//Calcolo con l'approssimazione di Gauss
  N = 4;
  int Ng = 2; //numero punti di Gauss
  
  double gauss = GaussRule(Root, c, d, 1, 3);
  cout << "Gauss:        " << scientific << setw(22) << setprecision(14) << gauss << ";" << endl;

//Adesso valuto la funzione Poly con Simpson e Gauss
  cout << "\nCalcolo una funzione diversa" << endl;

  cout << "\nExact:         " << -66./5. << ";" << endl;

//Calcolo con l'approssimazione di Simpson
  double e = -1.;
  double f = 5.;
  N = 2;

  double simpson = SimpsonRule(Poly, e, f, N);
  cout << "Simpson:      " << scientific << setw(22) << setprecision(14) << simpson << ";" << endl;
  
//Calcolo con l'approssimazione di Gauss
  gauss = GaussRule(Poly, e, f, 1, 3);
  cout << "Gauss:        " << scientific << setw(22) << setprecision(14) << gauss << ";" << endl;

  return 0;
}

double Func(double x){
  return exp(-x);
}

double Root(double x){
  return sqrt(1+x);
}

double Poly(double x){
  return 1 - x + 2*x*x + 0.5*x*x*x + 0.25*x*x*x*x - 0.125*x*x*x*x*x;
}

double RectangularRule(double (*F) (double), double a, double b, int N){

  double h = (b - a)/N; //larghezza delle divisioni dell'intervallo
  double sum = 0.; //dara' l'approssimazione dell'integrale
  double x; //ascissa di cui calcolo la funzione

  for(int i = 0; i < N; i++){
    x = a + i * h;
    sum += F(x);
//    cout << "X = " << x << endl;
  }

  sum *= h;

  return sum;
}

double TrapezioidalRule(double (*F) (double), double a, double b, int N){
  
  double h = (b - a) / N; //larghezza delle divisioni dell'intervallo
  double sum = 0.; //dara' l'approssimazione dell'integrale
  double x = a; //ascissa di cui calcolo la funzione
  double xp; //ascissa dell'iterazione precedente

  /*for(int i = 1; i <= N; i++){
    xp = x; 
    x = a + i * h;
    sum += 0.5 * (F(xp) + F(x));
//    cout << sum << endl;
  }*/


  sum = (F(a) + F(b))*0.5;

  for(int i=1; i<N; i++){
    x = a + i*h;
    sum += F(x);
  }
  
  sum *= h;

  return sum;
}

double SimpsonRule(double (*F) (double), double a, double b, int N){ //da determinare ancora
  
  double h = (b - a)/N;
  double sum = 0.;
  double x = a;  //ascissa in cui calcolo la funzione
  double w;  //pesi
  /*double f;
  
  for(int i=0;i<=N;i++){
    x = a + i*h;

    if(i==0||i==N)
      w = 1;
    else if(i%2==0 && i!=N && i!=0)
      w = 2;
    else 
      w = 4;
    f = (h*(F(x)*w))/3;
    sum += f;
  }

  return sum;*/

  sum = F(a) + F(b);
  w = 4.;
  
  for(int i=1; i<N; i++){
    x = a + i*h;
    sum += F(x)*w;
    w = 6. - w;
  }

  return sum*h/3.;
}
