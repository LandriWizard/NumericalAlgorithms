#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

void Product(double**, double*, double*, int);

int main(){
  cout << fixed << setprecision(4);

  int nrow = 4;
  int ncol = 4;
  int i,j;

  double **A;
  double *B = new double[nrow];
  double *C = new double[nrow];

  A    = new double*[nrow];
  A[0] = new double [ncol*nrow];
  for(i = 1; i < nrow; i++){
    A[i] = A[i - 1] + ncol;
  }
  A[0][0] = 1., A[0][1] =  3., A[0][2] = 2., A[0][3] = -4.;
  A[1][0] = 7., A[1][1] =  2., A[1][2] = 4., A[1][3] =  1.;
  A[2][0] = 0., A[2][1] = -1., A[2][2] = 2., A[2][3] =  2.;
  A[3][0] = 6., A[3][1] =  3., A[3][2] = 0., A[3][3] =  1.;

  for(i = 0; i < nrow; i++){
    for(j = 0; j < ncol; j++){
      cout << setw(10) << right << A[i][j] << " "; 
    }
    cout << endl;
  }

  cout << "---------------------" << endl;

  B[0] = 1., B[1] = 0., B[2] = 3., B[3] = 2.;
  for(i = 0; i < nrow; i++){
    cout << setw(10) << right << B[i] << endl; 
  }

  cout  << "--------------------" << endl;
 
  Product(A,B,C,4);
  for(i = 0; i < nrow; i++){
    cout << setw(10) << right << C[i] << endl; 
  }

  delete A[0];
  delete A;
  delete B;
  delete C;

  return 0;
}

void Product(double** A, double* X, double* Y, int n){

  int i, j;

  for(i = 0; i < n; i++){
    Y[i] = 0.;
    for(j = 0; j < n; j++)
      Y[i] += A[i][j]*X[j];
  }
}
