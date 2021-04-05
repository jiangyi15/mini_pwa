#include "tensor.h"
#include <cmath>
#include <iostream>

#ifndef DFUNCTION_H_
#define DFUNCTION_H_ 1

double factorial2(int x);

inline int min(int x, int y) {
  if (x > y) {
    return y;
  }
  return x;
}
inline int max(int x, int y) {
  if (x > y) {
    return x;
  }
  return y;
}

Tensor<double> small_d_weight(int j);

Tensor<double> small_d_function(int j, Tensor<double> *weight,
                                Tensor<double> *theta);

#endif
