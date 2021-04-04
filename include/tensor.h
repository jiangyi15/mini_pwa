#include <map>
#include <string>
#include <vector>

#ifndef M_PWATENSOR_H_
#define M_PWATENSOR_H_

typedef std::vector<int> Shape;

template <class T> class Tensor {
public:
  Tensor(Shape shape) : shape(shape) {
    int total_shape = 1;
    for (auto i : shape) {
      total_shape *= i;
    }
    this->ptr = (T *)malloc(total_shape * sizeof(T));
    for (int i = 0; i < total_shape; i++) {
      this->ptr[i] = 0;
    }
  }
  Tensor(int n) : shape({n}) {
    int total_shape = 1;
    for (auto i : shape) {
      total_shape *= i;
    }
    this->ptr = (T *)malloc(total_shape * sizeof(T));
    for (int i = 0; i < total_shape; i++) {
      this->ptr[i] = 0;
    }
  };
  ~Tensor() { free(this->ptr); }
  T *ptr;
  Shape shape;
};

typedef std::map<std::string, Tensor<double> *> MapTensor;

#endif
