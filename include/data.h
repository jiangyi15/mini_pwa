
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <complex>
#include <iostream>

#ifndef M_PWADATA_H_
#define M_PWADATA_H_

typedef std::vector<int> Shape;

template <class T> class Tensor {
public:
  Tensor(Shape shape) : shape(shape) {
    int total_shape = 1;
    for (auto i : shape) {
      total_shape *= i;
    }
    this->ptr = (T *)malloc(total_shape * sizeof(T));
  }
  ~Tensor() { free(this->ptr); }
  T *ptr;
  Shape shape;
};

typedef std::map<std::string, Tensor<double> *> MapTensor;

class EularAngle {
public:
  EularAngle(Tensor<double> *alpha, Tensor<double> *beta, Tensor<double> *gamma)
      : alpha(alpha), beta(beta), gamma(gamma){};
  EularAngle(std::vector<Tensor<double> *> data)
      : alpha(data[0]), beta(data[1]), gamma(data[2]){};
  Tensor<double> *alpha, *beta, *gamma;
};

class DecayData {
public:
  DecayData(EularAngle *angle) : angle(angle){};
  DecayData(EularAngle angle) : angle(&angle){};
  EularAngle *angle;
};

class ChainData {
public:
  ChainData(std::vector<DecayData *> data, MapTensor data_p)
      : data(data), data_p(data_p){};
  ChainData(std::pair<std::vector<DecayData *>, MapTensor> data)
      : data(data.first), data_p(data.second){};
  std::vector<DecayData *> data;
  std::map<std::string, Tensor<double> *> data_p;
};

#endif
