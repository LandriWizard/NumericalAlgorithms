#include "my_header.h"

#define KMAX 100 
#define ARRAYMAX 256

double Bisection(double (*F) (double), double a, double b, double xtol, double ytol, int &k, double &root){

  k = 0; //contatore delle iterazioni 
  double xm; //valore medio dell'intervallo
  double dx; //larghezza dell'intervallo
  double fa = F(a);
  double fm;

  do{
    xm = (a + b) / 2.;
    dx = fabs(a - b);
    k++;
    if(k > KMAX)
      break;
    fm = F(xm);
    //cout << "Bisection(): k = " << k << "; [a,b] = [" << a << ", " << b << "]; xm = " << xm << "; dx = " << dx << "; fm = " << fm << endl;
    if(fa*fm < 0)
      b = xm;
    else if(fa*fm > 0){
      a = xm;
      fa = fm;
    }
  }while(dx > xtol || fabs(fm) > ytol);
  root = xm;

  return 0;
}

double FalsePos(double (*F) (double), double a, double b, double tol, int &k, double &root){
  double fa, fb;
  double x0; //zero della retta
  double xp; //variabile in cui salvo l'x precedente
  double f0;
  double dx;
  k = 0;
  fa = F(a);
  fb = F(b);

  do{
    k++;
    if(k == KMAX) break;
    xp = x0;
    x0 = - ( b*fa - a*fb )/( fb - fa );
    f0 = F(x0);
    dx = fabs(xp - x0);
    //cout << "FalsePos(): k = " << k << "; [a,b] = [" << a << ", " << b << "]; x0 = " << x0 << "; |del| = " << dx << "; f0 = " << f0 << endl;
    if(fa*f0 < 0){
      b = x0;
      fb = f0;
    }
    else if(fa*f0 > 0){
      a = x0;
      fa = f0;
    }
  }while(dx > tol || fabs(f0) > tol);
  root = x0;

  return 0;
}

double Secant(double (*F) (double), double a, double b, double tol, int &k, double &root){
  double fa, fb;
  double dx;
  k = 0;
  fa = F(a);
  fb = F(b); 
  do{
    k++;
    if(k == KMAX) break;
    dx = fb*(b - a)/(fb - fa);
    a = b;
    fa = fb;
    b = b - dx;
    fb = F(b);
    //cout << "Secant(): k = " << k << "; xa = " << a << "; xb = " << b << "; dx = " << dx << endl;
  }while(fabs(dx) > tol);
  root = a;

  return 0;
}

double Newton(double (*F) (double), double (*G) (double), double a, double b, double tol, int &k,double &root){
  double xa = a, xb = b;
  double fa = F(xa), fb = F(xb);
  if(fa == 0.){
    k = 0;
    root = xa;
    return 0;
  }

  double xc, fc, dx;
  xc = 0.5*(xa + xb);
  if(fa*fb< 0.){
    for(int i = 1; i <128; i++){
      fc = F(xc);
      dx = fc/G(xc);
      xc -= dx;
      if(fabs(dx) < tol || fabs(fc) < tol){
        k = i;
        root = xc;
        return 0;
      }
    }
    return 1;
  }
  else return 1;

  return 0;
}

void Bracket(double (*F)(double), double a, double b, double *xL, double *xR, int n, int &nroots){
    double h = (b - a) / (double)n;
    double xi = a, xf = b , xnew;
    double fi = F(xi), ff;
    int k = 0;
    for(int i = 0; i < n; i++){
        ff = F(xi + (i + 1) * h);
        if(fi * ff < 0){
            xL[k] = xi + i * h;
            xR[k] = xi + (i + 1) * h;
            k++;
        }
        fi = ff;
    }
    nroots = k;
}
