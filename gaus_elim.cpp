#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

void SystemSolve(double**, double*, double*, int);
void SwapRows(double**, double*, int, int, int);

#define TRUE 1
#define FALSE 0
#define SESSION 3
#define DEBUG TRUE

int main(){
  cout << fixed << setprecision(4);

#if SESSION == 1
  int n = 3;
  int i;

  double **A;
  double *B = new double[n];
  double *X = new double[n];

  A    = new double*[n];
  A[0] = new double [n*n];
  for(i = 1; i < n; i++)
    A[i] = A[i - 1] + n;
 
  A[0][0] =  2., A[0][1] = -1., A[0][2] =  0.;
  A[1][0] = -1., A[1][1] =  2., A[1][2] = -1.;
  A[2][0] =  0., A[2][1] = -1., A[2][2] =  2.;

  B[0] = 1., B[1] = 2., B[2] = 1.;
 
  SystemSolve(A,B,X,n);
  for(i = 0; i < n; i++)
    cout << X[i] << endl;
  

  delete A[0];
  delete A;
  delete B;
  delete X;
#endif

#if SESSION == 2  
  int n = 4;
  int i;

  double **A;
  double *B = new double[n];
  double *X = new double[n];

  A    = new double*[n];
  A[0] = new double [n*n];
  for(i = 1; i < n; i++)
    A[i] = A[i - 1] + n;

  A[0][0] = 1., A[0][1] =  2., A[0][2] = 1., A[0][3] = -1.;
  A[1][0] = 3., A[1][1] =  2., A[1][2] = 4., A[1][3] =  4.;
  A[2][0] = 4., A[2][1] =  4., A[2][2] = 3., A[2][3] =  4.;
  A[3][0] = 2., A[3][1] =  0., A[3][2] = 1., A[3][3] =  5.;

  B[0] = 5., B[1] = 16., B[2] = 22., B[3] = 15.;

  SystemSolve(A,B,X,n);

  for(i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      cout << setw(10) << right << A[i][j] << " "; 
    }
    cout << endl;
  }

  cout << endl;

  for(i = 0; i < n; i++)
    cout << B[i] << endl;

  cout << endl;

  for(i = 0; i < n; i++)
    cout << X[i] << endl;
  
  delete A[0];
  delete A;
  delete B;
  delete X;
#endif

#if SESSION == 3 
  int n = 4;
  int i;

  double **A;
  double *B = new double[n];
  double *X = new double[n];

  A    = new double*[n];
  A[0] = new double [n*n];
  for(i = 1; i < n; i++)
    A[i] = A[i - 1] + n;

  A[0][0] = 1., A[0][1] =  2., A[0][2] = 1., A[0][3] = -1.;
  A[1][0] = 3., A[1][1] =  6., A[1][2] = 4., A[1][3] =  4.;
  A[2][0] = 4., A[2][1] =  4., A[2][2] = 3., A[2][3] =  4.;
  A[3][0] = 2., A[3][1] =  0., A[3][2] = 1., A[3][3] =  5.;

  B[0] = 5., B[1] = 16., B[2] = 22., B[3] = 15.;

  SystemSolve(A,B,X,n);

#if DEBUG == TRUE
  cout << "MATRICI A E B MODIFICATE" << endl;

  for(i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      cout << setw(10) << right << A[i][j] << " "; 
    }
    cout << endl;
  }

  cout << endl;

  for(i = 0; i < n; i++)
    cout << B[i] << endl;

  cout << endl;
#endif

  cout << "VETTORE SOLUZIONE" << endl;

  for(i = 0; i < n; i++)
    cout << X[i] << endl;
  
  delete A[0];
  delete A;
  delete B;
  delete X;
#endif


  return 0;
}

void SystemSolve(double **A, double *B, double *X, int n){
  int i,j,k,i_max;
  double g, temp, A_max, A_temp, B_temp;

  for(k = 0; k < n-1; k++){
//PARTIAL PIVOTING
    A_max = fabs(A[k][k]);
    i_max = k;
    for(i = k+1; i < n; i++){
      A_temp = fabs(A[i][k]);
      if(A_temp > A_max){
        i_max = i;
        A_max = A_temp;
      }
    }
/*    for(i = 0; i < n; i++){        //Row swappin
      temp        = A[k][i];
      A[k][i]     = A[i_max][i];
      A[i_max][i] = temp;
    }                              
    B_temp = B[k];
    B[k] = B[i_max];
    B[i_max] = B_temp;*/             //Rows swapped

    SwapRows(A, B, k, i_max, n);

//END PARTIAL PIVOTING

#if DEBUG == TRUE

    cout << "STAMPO OGNI MODIFICA ALLA MATRICE A E AL VETTORE B" << endl;

    for(i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        cout << setw(10) << right << A[i][j] << " "; 
      }
      cout << endl;
    }

    cout << endl;

    for(i = 0; i < n; i++)
      cout << B[i] << endl;

    cout << endl << endl;
#endif

    for(i = k+1; i < n; i++){
      g = A[i][k]/A[k][k];
      for(j = k+1; j < n; j++)
        A[i][j] -= g*A[k][j];
      A[i][k] = 0.;
      B[i] -= g*B[k];
    }
  }
 
  for(i = n-1; i >= 0; i--){
    temp = B[i];
    for(j = n-1; j>i; j--){
      temp -= X[j]*A[i][j];
    }
    X[i] = temp/A[i][i];
  }
}


void SwapRows(double** A, double* B, int k, int i_max, int n){
  double temp, B_temp;

  for(int i = 0; i < n; i++){        //Row swappin
    temp        = A[k][i];
    A[k][i]     = A[i_max][i];
    A[i_max][i] = temp;
  }                              
  B_temp = B[k];
  B[k] = B[i_max];
  B[i_max] = B_temp;             //Rows swapped
}
