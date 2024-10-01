#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#define NGMAX 32

double Func(double);
double Root(double);
double Poly(double);
double RectangularRule(double (*F) (double), double, double, int);
double TrapezioidalRule(double (*F) (double), double, double, int);
double SimpsonRule(double (*F) (double), double, double, int);
double GaussRule(double (*F) (double), double, double, int, int);

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

  /*while(err > tol){
    //double rectangular = RectangularRule(func, a, b, N);
    err = fabs(RectangularRule(Func, a, b, N) - RectangularRule(Func, a, b, N/2));
    N *= 2;
  }
    double rectangular = RectangularRule(Func, a, b, N);
    cout << "Rectangular:  " << scientific << setw(22) << setprecision(14) << rectangular << "; N = " << N/2 << endl;*/
  
  while(err > tol){
    sum0 = sum;
    sum = RectangularRule(Func, a, b, N);
    err = fabs(sum - sum0);
    N *= 2;
  }
    //double rectangular = RectangularRule(Func, a, b, N);
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
  //double trapezioidal = TrapezioidalRule(Func, a, b, N); 
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
  //double simpson = SimpsonRule(Func, a, b, N);
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
    //err = fabs(SimpsonRule(Root, c, d, N) - SimpsonRule(Root, c, d, N/2));
    err = fabs(sum - sum0);
    N *= 2;
  }
  //simpson = SimpsonRule(Root, c, d, N);
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
 
  gauss = GaussRule(Poly, e, f, 1, 4);
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

double GaussRule(double (*F)(double), double a, double b, int N, int Ng){
  
  double sum = 0.;
  double sumk;
  double x[NGMAX]; //array con le ascisse
  double w[NGMAX]; //array con i pesi
  double xv1;
  double xv2;
  double wv1;
  double wv2;
  double h = (b - a)/N;
  double x0;
  double x1;
  double deltax;
  double avx;
  
//riempio gli array in base al numero di punti di Gauss
  if(Ng == 2){
    xv1 = 1./sqrt(3.);
    x[0] = -xv1; x[1] = xv1;
    w[0] = 1; w[1] = 1;
  }
  else if(Ng == 3){
    xv1 = sqrt(3./5.);
    x[0] = -xv1; x[1] = 0; x[2] = xv1;
    w[0] = 5./9.; w[1] = 8./9.; w[2] = 5./9.;
  }
  else if(Ng == 4){
    xv1 = 0.339981043584856;
    xv2 = 0.861136311594053;
    x[0] = -xv1; x[1] = xv1; x[2] = -xv2; x[3] = xv2;
    wv1 = 0.652145154862546;
    wv2 = 0.347854845137454;
    w[0] = wv1; w[1] = wv1; w[2] = wv2; w[3] = wv2;
  }

  for(int i = 0; i < N; i++){
    x0 = a + i*h;
    x1 = a + (i+1)*h;
    deltax = 0.5*(x1 - x0);
    avx = 0.5*(x1 + x0);
    sumk = 0.;
    for(int k = 0; k < Ng; k++){
      sumk += w[k]*F(deltax*x[k] + avx);
    }
    sumk *= deltax;
    sum += sumk;
  }

  return sum;
}
