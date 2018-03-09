#include <iostream>
using namespace std;

class proof{
  int a;
  int b;
public:
  proof(int);
  void example();
};

proof::proof(int c){
  a = 0;
  a = 1 + c;
}

void proof::example(){
  cout << "a: " << a;
}

int main(){
  int l = 5;
  proof p(l);
  p.example();
  return 0;
}
