#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
using namespace std;

int main(){
  
  cout << setiosflags(ios::scientific);
  cout << setprecision(14);

  double N = 1.e3; //numero di atomi, varia con il tempo
  double Ni = N;   //numero iniziale di atomi
  double M;        //numero di atomi normalizzato al numero iniziale
  double Nd;       //numero di atomi decaduti per ciclo
  int dt = 0; 
  int tmax = 500;  //tempo massimo
  double lambda = 0.01; //probabilità di decadimento
  double r;//numero randomico

  srand48(time(NULL));  //seed del rand


//definisco il file che manderò a gnuplot
  ofstream fdata;
  fdata.open("decay.dat");

  while(N > 0){

    if(dt == tmax)  //definisco un limite di massimo di tempo
      break;

    Nd = 0; //inizializzo il numero di atomi decaduti in questo ciclo a zero

    for(int i=0; i<N; i++){  //per ogni atomo calcolo la probabilità di decadimento

      r = drand48();

      if(r<lambda)  //conto gli atomi decaduti
        Nd++;
    }

    N -= Nd;     //tolgo gli atomi decaduti in questo ciclo
    M = N/Ni;    //normalizzo il numero di atomi
    fdata << dt << " " << M << endl;  //riempio il file che manderò a gnuplot
    dt++;  //aumento il counter del tempo
  }

  fdata.close();

  return 0;
}
