#include <iostream>
using namespace std;

int main(){
  int guess = 0;
  int n;
  int li = 1; //limite inferiore
  int ls = 100; //limite superiore
  int c = 1; //contatore di guess fatte

  srand(time(NULL));
  n = rand()%100 + 1;

  while(guess != n){
    cout << "n in [" << li << ", " << ls << "]" << endl;
    cout << "type your guess #" << c << " > ";
    cin >> guess; 
    cout << endl;
    if(guess < n){
      li = guess;
      c++;
    }
    else if(guess > n){
      ls = guess;
      c++;
    }
  }
  
  cout << "Number found in " << c << " iterations !!" << endl;

  return 0;
}
