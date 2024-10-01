#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

#define IBEG 1
#define IEND NX-2
#define JBEG 1
#define JEND NY-2

#define SESSION 3

#define JACOBI 1
#define GAUSS_SIEDEL 2
#define MSOR 3
#define METHOD MSOR


double SourceTerm(double, double);
double ExactValue(double, double);
void BoundaryValues(double**, double*, double*, int, int);
void Jacobi(double**, double**, double*, double*, double, double, int, int, int&, double&);
void GaussSiedel(double**, double*, double*, double, double, int, int, int&, double&);
void SOR(double**, double*, double*, double, double, int, int, int&, double&);

double g_a = 0.1;
double g_rho = 1.;

int main(){
  cout << setiosflags(ios::scientific);
  cout << setprecision(8);
  cout << setw(2);
  
#if SESSION == 1
  int i,j,k = 0;
  int NX = 32;
  int NY = 32; 
  double xi = 0., xf = 1.;
  double yi = 0., yf = 1.;
  double h = fabs(xf - xi)/(double)(NX - 1); //same value for y
  double err;
  double tol = 1.e-7;

  double **F0;
  double **F1;
  double *X = new double[NX];
  double *Y = new double[NY];

  F0    = new double*[NY];
  F0[0] = new double [NX*NY];
  for(i = 1; i < NY; i++){
    F0[i] = F0[i - 1] + NX;
  }

  F1    = new double*[NY];
  F1[0] = new double [NX*NY];
  for(i = 1; i < NY; i++){
    F1[i] = F1[i - 1] + NX;
  }

//setting values in coordinates arrays
  for(i = 0; i < NX; i++) X[i] = xi + i*h;
  for(i = 0; i < NY; i++) Y[i] = yi + i*h;

#if METHOD == JACOBI
// Jacobi cycle
  k = 0;
  err = 0.;

  Jacobi(F0,F1,X,Y,h,tol,NX,NY,k,err);

  cout << "JACOBI" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("jacobi.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F1[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif

#if METHOD == GAUSS_SIEDEL
//Gauss Siedel cycle
  k = 0;
  err = 0.;

  GaussSiedel(F0,X,Y,h,tol,NX,NY,k,err);

  cout << "GAUSS-SIEDEL" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("gauss_siedel.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F0[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif

#if METHOD == MSOR
//SOR cycle
  k = 0;
  err = 0.;

  SOR(F0,X,Y,h,tol,NX,NY,k,err);

  cout << "SOR" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("sor.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F0[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif

#if SESSION == 2
  int i,j,k;
  int NX = 128;
  int NY = 128; 
  double xi = -1., xf = 1.;
  double yi = -1., yf = 1.;
  double h = fabs(xf - xi)/(double)(NX - 1); //same value for y
  double err;
  double tol = 1.e-7;

  double **F0;
  double **F1;
  double *X = new double[NX];
  double *Y = new double[NY];

  F0    = new double*[NY];
  F0[0] = new double [NX*NY];
  for(i = 1; i < NY; i++){
    F0[i] = F0[i - 1] + NX;
  }

  F1    = new double*[NY];
  F1[0] = new double [NX*NY];
  for(i = 1; i < NY; i++){
    F1[i] = F1[i - 1] + NX;
  }
#endif

//setting values in coordinates arrays
  for(i = 0; i < NX; i++) X[i] = xi + i*h;
  for(i = 0; i < NY; i++) Y[i] = yi + i*h;

#if METHOD == JACOBI
// Jacobi cycle
  k = 0;
  err = 0.;

  Jacobi(F0,F1,X,Y,h,tol,NX,NY,k,err);

  cout << "JACOBI" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("jacobi.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F1[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif

#if METHOD == GAUSS_SIEDEL
//Gauss Siedel cycle
  k = 0;
  err = 0.;

  GaussSiedel(F0,X,Y,h,tol,NX,NY,k,err);

  cout << "GAUSS-SIEDEL" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("gauss_siedel.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F0[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif

#if METHOD == MSOR
//SOR cycle
  k = 0;
  err = 0.;

  SOR(F0,X,Y,h,tol,NX,NY,k,err);

  cout << "SOR" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("sor.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F0[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif
#endif

#if SESSION == 3
  int i,j,k;
  int NX = 129;
  int NY = 65; 
  double xi = 0., xf = 2.;
  double yi = 0., yf = 1.;
  double h = fabs(xf - xi)/(double)(NX - 1); //same value for y
  double err;
  double tol = 1.e-7;

  double **F0;
  double **F1;
  double *X = new double[NX];
  double *Y = new double[NY];

  F0    = new double*[NY];
  F0[0] = new double [NX*NY];
  for(i = 1; i < NY; i++){
    F0[i] = F0[i - 1] + NX;
  }

  F1    = new double*[NY];
  F1[0] = new double [NX*NY];
  for(i = 1; i < NY; i++){
    F1[i] = F1[i - 1] + NX;
  }

//setting values in coordinates arrays
  for(i = 0; i < NX; i++) X[i] = xi + i*h;
  for(i = 0; i < NY; i++) Y[i] = yi + i*h;

#if METHOD == JACOBI
// Jacobi cycle
  k = 0;
  err = 0.;

  Jacobi(F0,F1,X,Y,h,tol,NX,NY,k,err);

  cout << "JACOBI" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("jacobi.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F1[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif

#if METHOD == GAUSS_SIEDEL
//Gauss Siedel cycle
  k = 0;
  err = 0.;

  GaussSiedel(F0,X,Y,h,tol,NX,NY,k,err);

  cout << "GAUSS-SIEDEL" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("gauss_siedel.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F0[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif

#if METHOD == MSOR
//SOR cycle
  k = 0;
  err = 0.;

  SOR(F0,X,Y,h,tol,NX,NY,k,err);

  cout << "SOR" << endl;
  cout << "#iterazioni = " << k-1 << "; errore = " << err << endl;

  ofstream fdata;
  fdata.open("sor.dat");
  for(i = 0; i < NX; i++){
    for(j = 0; j < NY; j++){
      fdata << X[i] << " " << Y[j] << " " << F0[i][j] << endl;
    }
    fdata << endl;
  }
  fdata.close();
#endif
#endif

  delete F0[0];
  delete F0;  
  delete F1[0];
  delete F1;  
  delete X;
  delete Y;
  return 0;
}

#if SESSION == 1
double SourceTerm(double x, double y){
  return 0.;
}

double ExactValue(double x, double y){
  double S = SourceTerm(x,y);
  return exp(-M_PI*x)*sin(-M_PI*y) + 0.25*S*(x*x + y*y);
}
#endif

#if SESSION == 2
double SourceTerm(double x, double y){
  double r = sqrt(x*x + y*y);

  if(r <= g_a)
    return -g_rho;
  else 
    return 0.;
}

double ExactValue(double x, double y){
  double r = sqrt(x*x + y*y);
  double S = SourceTerm(x,y);

  if(r <= g_a)
    return -0.25*g_rho*r*r;
  else
    return -0.5*g_rho*g_a*g_a*(log(r/g_a) + 0.5);
}
#endif

#if SESSION == 3
double SourceTerm(double x, double y){
  return 0.;
}
#endif

void BoundaryValues(double **M, double *X, double *Y, int NX, int NY){
#if SESSION == 1 || SESSION == 2
  for(int i = 0; i < NX; i++){
    M[i][0]    = ExactValue(X[i],Y[0]);    
    M[i][NY-1] = ExactValue(X[i],Y[NY-1]);    
  }
  for(int j = 0; j < NY; j++){
    M[0][j]    = ExactValue(X[0],Y[j]);    
    M[NX-1][j] = ExactValue(X[NX-1],Y[j]);    
  }
#endif

#if SESSION == 3
  double sigma, h;
  int i, j;
  h = (X[NX - 1] - X[0])/(double)(NX - 1);

//bottom: Dirichlet
  j = 0;
  for(i = 0; i < NX; i++){
    M[i][j]    = 0.;    
  }

//top: Dirichlet
  j = NY - 1;
  for(i = 0; i < NX; i++){
    M[i][j] = 2. - X[i];    
  }

//left: Neumann
  sigma = 0.;
  i = 0;
  for(j = 0; i < NX; i++){
    M[i][j] = M[i+1][j] - h*sigma;    
  }

//right: Neumann
  sigma = 3.;
  i = NX - 1;
  for(j = 0; i < NX; i++){
    M[i][j] = M[i-1][j] + h*sigma;    
  }  
#endif
}

void Jacobi(double **F0, double **F1, double *X, double *Y, double h, double tol, int NX, int NY, int &k, double &err){

  int i, j;
  double deltax, deltay;

  for(i = IBEG; i <= IEND; i++){
    for(j = JBEG; j <= JEND; j++){
      F0[i][j] = 0.;
      F1[i][j] = 0.;
    }
  }

//Jacobi cycle
  do{
    k++;
    err = 0.;

    BoundaryValues(F0,X,Y,NX,NY);
    BoundaryValues(F1,X,Y,NX,NY);

//computing new internal values
    for(i = IBEG; i <= IEND; i++){
      for(j = JBEG; j <= JEND; j++){
        F1[i][j] = 0.25*(F0[i+1][j] + F0[i-1][j] + F0[i][j+1] + F0[i][j-1] - h*h*SourceTerm(X[i],Y[j]));
      } 
    }

//computing residual
    for(i = IBEG; i <= IEND; i++){
      for(j = JBEG; j <= JEND; j++){
        deltax = F1[i+1][j] - 2.*F1[i][j] + F1[i-1][j];
        deltay = F1[i][j+1] - 2.*F1[i][j] + F1[i][j-1];
        err += fabs(deltax + deltay - h*h*SourceTerm(X[i],Y[j]));
      } 
    }

//update matrix
    for(i = IBEG; i <= IEND; i++){
      for(j = JBEG; j <= JEND; j++){
        F0[i][j] = F1[i][j];
      } 
    }

    cout << endl;
    cout << k << "\t" << err << endl;

  }while(err > tol);
}

void GaussSiedel(double **F0, double *X, double *Y, double h, double tol, int NX, int NY, int &k, double &err){
  int i, j;
  double deltax, deltay;

  for(i = IBEG; i <= IEND; i++){
    for(j = JBEG; j <= JEND; j++){
      F0[i][j] = 0.;
    }
  }

  do{
    k++;
    err = 0.;

    BoundaryValues(F0,X,Y,NX,NY);

    for(i = IBEG; i <= IEND; i++){
      for(j = JBEG; j <= JEND; j++){
        F0[i][j] = 0.25*(F0[i+1][j] + F0[i-1][j] + F0[i][j+1] + F0[i][j-1] - h*h*SourceTerm(X[i],Y[j]));
      } 
    }

    for(i = IBEG; i <= IEND; i++){
      for(j = JBEG; j <= JEND; j++){
        deltax = F0[i+1][j] - 2.*F0[i][j] + F0[i-1][j];
        deltay = F0[i][j+1] - 2.*F0[i][j] + F0[i][j-1];
        err += fabs(deltax + deltay - h*h*SourceTerm(X[i],Y[j]));
      } 
    }

    cout << endl;
    cout << k << "\t" << err << endl;

  }while(err > tol);
}

void SOR(double **F0, double *X, double *Y, double h, double tol, int NX, int NY, int &k, double &err){
  int i, j;
  double deltax, deltay;
  double omega = 2./(1. + (M_PI/(double)NX));

//  cout << "SOR\n\n";

  for(i = IBEG; i <= IEND; i++){
    for(j = JBEG; j <= JEND; j++){
      F0[i][j] = 0.;
    }
  }

  do{
    k++;
    err = 0.;

    BoundaryValues(F0,X,Y,NX,NY);

    for(i = IBEG; i <= IEND; i++){
      for(j = JBEG; j <= JEND; j++){
        F0[i][j] = (1. - omega)*F0[i][j] + 0.25*omega*(F0[i-1][j] + F0[i+1][j] + F0[i][j-1] + F0[i][j+1] - h*h*SourceTerm(X[i],Y[j]));
      } 
    }

    for(i = IBEG; i <= IEND; i++){
      for(j = JBEG; j <= JEND; j++){
        deltax = F0[i+1][j] - 2.*F0[i][j] + F0[i-1][j];
        deltay = F0[i][j+1] - 2.*F0[i][j] + F0[i][j-1];
        err += fabs(deltax + deltay - h*h*SourceTerm(X[i],Y[j]));
      } 
    }

    cout << endl;
    cout << k << "\t" << err << endl;

  }while(err > tol);
}
