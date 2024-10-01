#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

void TridSolver(double*,double*,double*,double*,double*,int);

int main(){
  int i;
  double x;          //coordinata x
  double xi = 0.;    //x iniziale
  double xf = 10.;   //x finale
  int n = 100;       //numero di punti della griglia
  double dx = fabs(xf - xi)/(double)(n-1);  //distanza tra due punti
  double YL = 0.;    //Y(xi)
  double YR = 1.;    //Y(xf)

  double *A;
  double *B;
  double *C;
  double *R;
  double *Y;
  A = new double[n];
  B = new double[n];
  C = new double[n];
  R = new double[n];
  Y = new double[n];
	
  Y[0] = YL;     //assegno le condizioni a contorno a Y
  Y[n-1] = YR;
	
  for(i = 0; i < n; i++){
    x = xi + i*dx;                    //aggiornamento della posizione x
    A[i] = x*(x - 0.5*dx);            //coefficiente di Y[i-1]
    B[i] = x*x*(dx*dx - 2.) - dx*dx;  //coefficiente di Y[i]
    C[i] = x*(x + 0.5*dx);            //coefficiente di Y[i+1]
    R[i] = 0.;                        //risultato della i-esima equazione
  }

  R[1] = R[1] - A[1]*Y[0];            
  R[n-2] = R[n-2] - C[n-2]*Y[n-1];    
	
  TridSolver(A+1,B+1,C+1,R+1,Y+1,n-2);   //risolve il sistema per Y
	
  ofstream fdata;
  fdata.open("bessel.dat");	
  for(i = 0; i < n; i++){
    x = xi + i*dx;
    fdata << x << " " << Y[i] << endl;
  }
  fdata.close();
  
  delete A;	
  delete B;	
  delete C;	
  delete R;	
  delete Y;	
  return 0;
}

void TridSolver(double* A, double* B, double* C, double* R, double* X, int n){
  int i;
  double *H;
  double *P;
  H = new double[n];
  P = new double[n];
  H[0] = C[0]/B[0];
  P[0] = R[0]/B[0];
  for(i = 1; i < n; i++){
    H[i] = C[i]/(B[i] - A[i]*H[i-1]);
    P[i] = (R[i] - A[i]*P[i-1])/(B[i] - A[i]*H[i-1]);
  }
  X[n-1] = P[n-1];
  for(i = n-2; i >= 0; i--){
    X[i] = P[i] - H[i]*X[i+1];
  }
  delete P;
  delete H;
}
