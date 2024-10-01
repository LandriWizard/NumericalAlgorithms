#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

void SwapRows(double**,double*,int,int,int);
void PrintVector(double*,int);
void PrintMatrix(double**,int);
void TridSolver(double*,double*,double*,double*,double*,int);

int main(){
  cout << setiosflags(ios::scientific);
  cout << fixed << setprecision(4); 

  int n = 32;  

  double *A;
  A = new double[n];
  double *B;
  B = new double[n];
  double *C;
  C = new double[n];
  double *R;
  R = new double[n];
  double *X;
  X = new double[n];
  double *Y;
  Y = new double[n];
  double *P;
  P = new double[n];
  double *H;
  H = new double[n];

  Y[0]=0.0;
  Y[n-1]=0.0;
  X[0]=0.0;
  X[n-1]=1.0;
  A[0]=0.0;
  C[n-1]=0.0;
  double dx = 1./double(n-1);

  for(int i=1; i<n; i++)    A[i] =    1.0;
  for(int i=0; i<n; i++)    B[i] =   -2.0;
  for(int i=0; i<n-1; i++)  C[i] =    1.0;
  for(int i=0; i<n; i++)    X[i] =   dx*i;
  for(int i=0; i<n; i++)    R[i] =  dx*dx;

  R[0] = 0;
  R[n - 1] = 0;
  
  TridSolver(A+1, B+1, C+1, R+1, Y+1, n-2);
 
  ofstream fdata;
  fdata << setiosflags(ios::scientific);
  fdata << fixed << setprecision(4); 
  fdata.open("bvp.dat");
  for(int i = 0; i < n; i++){
    fdata << X[i] << " " << Y[i] << endl;
  }
  fdata.close();
	
  // printing the result
//  PrintVector(Y,n);

  // deleting the allocated memory 
  delete A;
  delete B;
  delete C;
  delete R;
  delete X;
  delete Y;
  delete P;
  delete H;
	
  return 0;
}

void SwapRows(double **A, double *B, int r1, int r2, int n){
  double temp;
  // il seguente scambio si può riscrivere con "memcpy"
  for(int j=0;j<n;j++){
    temp = A[r1][j];
    A[r1][j] = A[r2][j];
    A[r2][j] = temp;
  }
  temp = B[r1];
  B[r1] = B[r2];
  B[r2] = temp;
}

void PrintVector(double *A, int n){
  for (int i = 0; i < n; i++){
    cout << setw(7) << right << A[i] << endl;
  }
}

void PrintMatrix(double **A, int n){
  cout << fixed << setprecision(4);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      cout << setw(10) << right << A[i][j] << " ";
    }
    cout << endl;
  }
}

void TridSolver(double *A, double *B, double *C, double *R, double *X, int n){
  double *H;
  H = new double[n];
  double *P;
  P = new double[n];
  H[0] = C[0]/B[0];
  P[0] = R[0]/B[0];
  for(int i = 1; i < n; i++){
    H[i] = C[i]/(B[i] - A[i]*H[i-1]);
    P[i] = (R[i] - A[i]*P[i-1])/(B[i] - A[i]*H[i-1]);
  }
  X[n-1] = P[n-1];
  for(int i = n-2; i >= 0; i--){
    X[i] = P[i] - H[i]*X[i+1];
  }
  delete P;
  delete H;
}
