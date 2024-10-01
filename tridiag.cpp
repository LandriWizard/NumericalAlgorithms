#include <iostream>
#include <iomanip>
#include <cmath>

void SwapRows(double**,double*,int,int,int);
void PrintVector(double*,int);
void PrintMatrix(double**,int);
void TridSolver(double*,double*,double*,double*,double*,int);

int main(){
  using namespace std;

  int N = 5; // order of the matrix

  double *A;
  A = new double[N];
  double *B;
  B = new double[N];
  double *C;
  C = new double[N];
  double *R;
  R = new double[N];
  double *X;
  X = new double[N];

  A[1] = 1;
  A[2] = 1;
  A[3] = 1;
  A[4] = 1;

  B[0] = 2;
  B[1] = 2;
  B[2] = 2;
  B[3] = 2;
  B[4] = 2;

  C[0] = 1;
  C[1] = 1;
  C[2] = 1;
  C[3] = 1;

  R[0] = 1;
  R[1] = 0;
  R[2] = 3;
  R[3] = 1;
  R[4] = 0;

  TridSolver(A,B,C,R,X,N);
	
  // printing the result
  PrintVector(X,N);

  // deleting the allocated memory 
  delete A;
  delete B;
  delete C;
  delete R;
  delete X;
	
  return 0;
}

void SwapRows(double **A, double *b ,int r1, int r2, int n){
  double tmp;
  // il seguente scambio si può riscrivere con "memcpy"
  for(int j=0;j<n;j++){
    tmp = A[r1][j];
    A[r1][j] = A[r2][j];
    A[r2][j] = tmp;
  }
  tmp = b[r1];
  b[r1] = b[r2];
  b[r2] = tmp;
}

void PrintVector(double *a, int N){
  using namespace std;
  for (int i = 0; i < N; i++){
    cout << setw(7) << right << a[i] << endl;
  }
}

void PrintMatrix(double **A, int N){
  using namespace std;
  cout << fixed << setprecision(4);
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      cout << setw(10) << right << A[i][j] << " ";
    }
    cout << endl;
  }
}

void TridSolver(double *A, double *B, double *C, double *R, double *X, int N){
  double *H;
  H = new double[N];
  double *P;
  P = new double[N];
  H[0] = C[0]/B[0];
  P[0] = R[0]/B[0];
  for(int i = 1; i<N; i++){
    H[i] = C[i]/(B[i]-A[i]*H[i-1]);
    P[i] = (R[i]-A[i]*P[i-1])/(B[i]-A[i]*H[i-1]);
  }
  X[N-1] = P[N-1];
  for(int i = N-2; i>=0; i--){
    X[i] = P[i]-H[i]*X[i+1];
  }
  delete P;
  delete H;
}
