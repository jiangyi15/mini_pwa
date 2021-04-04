
#include "tensor.h"
#include <complex>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>
#ifndef M_PWADATA_H_
#define M_PWADATA_H_

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
