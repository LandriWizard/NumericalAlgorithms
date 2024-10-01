#include <iostream>
#include <cmath>

void Average (double *x, int size, double&, double&, double&);  

int main()
{
  using namespace std;
  int i, n = 1000;
  double x[n], xav, s2, s;

  srand48(time(NULL));  // creo seed per sequenza di numeri casuali
 
  for (i = 0; i < n; i++) {
    
    x[i] = drand48();  // riempio con numero casuale tra 0 e 1

  }

  Average (x, n, xav, s2, s);
  cout << "<x> = " << xav << "; s^2 = " << s2 << "; s = " << s << endl;
  
  return 0;
  
}
// * puntatore è uguale a scrivere [] qui
void Average (double *x, int size, double& mu, double& s2, double& s)

{
  int i;
  
// Compute mean

  mu=0.0;
  for (i=0; i<size; i++){
    mu += x[i];
  }
  mu /= (double)size;

// Compute Variance
//due cicli for per forza perchè mi serve mu totale

  s2 = 0.0;
  for (i = 0; i < size; i++){
    s2 += (x[i] - mu)*(x[i] - mu);
  }
  s2 /= (double)size;

// Compute Standard deviation 

  s = sqrt(s2);

}