#include <iostream>

int main(){

  //Method #1, for

  std::cout << "Using for\n" << std::endl;

  for(int i=1;i<=10;i++){ //ciclo for che conta da 1 a 10          posso dichiarare la variabile direttamente nelle () del for
    if(i%2 == 1) std::cout << i << std::endl;
  }
  std::cout << "\nEnding for\n" << std::endl;

  //Method #2, while
 
  std::cout << "Using while\n" << std::endl;
  
  int j=1; //inizializzo qui la variabile usata in while perche non la posso dichiarare nelle () del while!!

  while(j<=10){ //ciclo while che conta da 1 a 10
    if(j%2 == 1) std::cout << j << std::endl;
    j++;
  }
  std::cout << "\nEnding while\n" << std::endl;

  return 0;

}
