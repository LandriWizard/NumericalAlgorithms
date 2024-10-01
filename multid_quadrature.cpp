#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#define NGMAX 32

static double g_ycoord; //valore fissato di y su cui integriamo le x
static int    g_nint = 1;   //intervalli della funzione quadratura

double Func2D(double, double); //la funzione 2D
double Func2DCircle(double); //altra funzione 2D
double FFunc1D(double); //la uso per integrare 1D in y
double FFunc1DCircle(double);
double GFunc1D(double); //la uso per integrare in 1D su x
double GFunc1DCircle(double);
double GaussRule(double (*F) (double), double, double, int, int);

int main(){

  double sum = 0.;
  double ybeg = -1., yend = 1.; //bordi del dominio in y
  
  sum = GaussRule(GFunc1D, ybeg, yend, g_nint, 3);

  cout << "Exact :" << scientific << setw(22) << setprecision(14) << 412./45. << ";" << endl;
  cout << "Gauss :" << scientific << setw(22) << setprecision(14) << sum      << ";" << endl;

//adesso voglio integrare la funzione Func2DCircle
  
  cout << "Inegro un cerchio unitario" << endl;

  sum = GaussRule(GFunc1DCircle, ybeg, yend, g_nint, 4);
  cout << "\nGauss :" << scientific << setw(22) << setprecision(14) << sum << ";" << endl;

  return 0;
}

double Func2D(double x, double y){

  return x*x*x*x*y*y + 2*x*x*y*y - x*x*y + 2;
}

double FFunc1D(double x){

  return Func2D(x, g_ycoord);
}

double GFunc1D(double y){
  double sum_x, xbeg = -1., xend = 1.;
  
  g_ycoord = y; 

  sum_x = GaussRule(FFunc1D, xbeg, xend, g_nint, 3);
  return sum_x;
}

double Func2DCircle(double x, double y){
  double r = sqrt(x*x + y*y);
  if(r <=1) 
    return 1;
  else
    return 0;
}

double FFunc1DCircle(double x){

  return Func2DCircle(x, g_ycoord);
}

double GFunc1DCircle(double y){
  double sum_x, xbeg = -1., xend = 1.;
  
  g_ycoord = y; 

  sum_x = GaussRule(FFunc1DCircle, xbeg, xend, g_nint, 3);
  return sum_x;
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
