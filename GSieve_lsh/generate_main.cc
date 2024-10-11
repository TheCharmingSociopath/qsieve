#include <iostream>
#include "NTL/LLL.h"

NTL_CLIENT

int main(int argc, char** argv) {
  long i;
  long dims = atoi(argv[1]);
  long temp = time(NULL);
  srand(temp);
  SetSeed(to_ZZ(temp));
  ZZ max;
  GenPrime(max, 10 * dims);
  mat_ZZ B;
  B.SetDims(dims,dims);
  for(i = 1;i < dims; i++){
    B[i][i] = 1;
    B[i][0] = RandomBnd(max);
  }
  B[0][0] = max;
  cout << B;
}