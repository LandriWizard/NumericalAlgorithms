#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#define nMAX 256

using namespace std;

int Newton(double (*func)(double), double (*dfunc)(double), double, double, double, int &, double &);
int Bisection(double (*func)(double), double, double, double, int &, double &);
void Bracketing(double (*func)(double), double, double, double *, double *, int , int &);


double func(double);
double dfunc(double);

double func2(double);
double dfunc2(double);

double func5(double);
double dfunc5(double);



int main(){
    /*
    double a = -5.0, b = 0.0, xtol = 1.e-8, root = 0.0;
    int ntry = 0;
    int r = Newton(func, dfunc, a, b, xtol, ntry, root);
    cout << "---- root = " << root << ", in " << ntry << " steps" << endl << endl;
    ntry = 0, root = 0.0;
    r = Newton(func2, dfunc2, a, b, xtol, ntry, root);
    cout << "---- root = " << root << ", in " << ntry << " steps" << endl; 
    */

    double a = -10.0, b = 10.0, xtol = 1.e-8, root = 0.0;
    int n = 10, nroots, r, ntry;
    double xL[nMAX], xR[nMAX];

    Bracketing(func5, a, b, xL, xR, n, nroots);

    for(int i = 0; i < nroots; i++){
        cout << "#" << i+1 << "; xL = " << xL[i] << "; xR = " << xR[i]<< endl;
    }
 
    for(int i = 0; i < nroots; i++){
        int r = Newton(func5, dfunc5, xL[i], xR[i], xtol, ntry, root);
        cout << root << endl;
    }
}


double func5(double x){
    return sin(x) - ((x*x/100.0) + x/5.0 + 1.0/3.0);
}
 double dfunc5(double x){
     return cos(x) - (x/50.0 + 1.0/5.0);
 }

double func(double x){
    return x*x*x - 3 * x*x + x + 5;
}

double dfunc(double x){
    return 3 * x*x - 6 * x + 1;
}

double func2(double x){
    int n = 3;
    double a[] = {5.0, 1.0, -3.0, 1.0};
    double p = a[n];
    double dp = 0.0;
    for(int i = n-1; i >= 0; i--){
        p = p * x + a[i];
    }
    return p;
}

double dfunc2(double x){
    int n = 3;
    double a[] = {5.0, 1.0, -3.0, 1.0};
    double p = a[n];
    double dp = 0.0;
    for(int i = n-1; i >= 0; i--){
        dp = dp * x + p;
        p = p * x + a[i];
    }
    return dp;
}


int Newton(double (*func)(double), double (*dfunc)(double) , double a0, double b0, double xtol, int &ntry, double &root){
    double xa = a0, xb = b0;
    double fa = func(xa), fb = func(xb);
    if(fa == 0.0){
        ntry = 0;
        root = xa;
        return 0;
    }
    if(fb == 0.0){
        ntry = 0;
        root = xb;
        return 0;
    }

    double xc, fc, dx;
    xc = 0.5*(xa + xb);
    if(fa * fb < 0.0){
        for(int i = 1; i < 128; i++){
            fc = func(xc);
            dx = fc / dfunc(xc);
            xc -= dx;
            // cout << "Newton(); i = " << i <<  "; xc = " << xc << " dx = "  << dx << " fc = " << fc << endl;
            if (fabs(dx) < xtol){
                ntry = i;
                root = xc;
                return 0;
            }
        }
        return 1;
    }else return 1;
}

int Bisection(double (*func)(double), double a, double b, double tol, int &ntry, double &zero){
    double xin = a;
    double xfin = b;
    double fi = func(a);
    double ff = func(b);

    if(fi * ff>0) return 1;

    for(int k=1; 1; k++){
        double xm = (xin + xfin)*0.5;
        double dx = fabs(xfin - xin);
        double fm = func(xm);
        if (dx < tol) break;
        if( fi * fm < 0. ){
            xfin = xm;
            ff = fm;
        }
        else if ( fi * fm > 0. ){ 
            xin = xm;
            fi = fm;
        }

        cout << " Bisection : " <<  " k = " << k << ";    " << " [a,b] = " << " [ " << xin << " , " << xfin << " ] " << " xm = " << xm << "    dx = "  << dx << "     fm = " << endl;

    }
    zero = xin;
    return 0;
}


void Bracketing(double (*func)(double), double a, double b, double *xL, double *xR, int n, int &nroots){
    double eps = (b - a) / (double)n;
    double xin = a, xfin = b , xnew;
    double fi = func(xin), ff;
    int k = 0;
    for(int i = 0; i < n; i++){
        ff = func(xin + (i + 1) * eps);
        if(fi * ff < 0){
            xL[k] = xin + i * eps;
            xR[k] = xin + (i + 1) * eps;
            k++;
        }
        fi = ff;
    }
    nroots = k;
}






