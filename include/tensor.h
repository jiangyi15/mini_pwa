#include <iostream>
#include <map>
#include <string>
#include <vector>
#ifndef M_PWATENSOR_H_
#define M_PWATENSOR_H_

typedef std::vector<int> Shape;

template <class T> class Tensor {
  bool _has_init = true;

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
  Tensor(Shape shape, T *ptr) : shape(shape), ptr(ptr), _has_init(false) {}
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
  Tensor<T> operator[](int n) {
    Shape shape;
    for (int i = 1; i < this->shape.size(); i++) {
      shape.push_back(this->shape[i]);
    }
    return Tensor<T>(shape,
                     this->ptr + n * this->total_shape() / this->shape[0]);
  };
  int total_shape() {
    int total_shape = 1;
    for (auto i : shape) {
      total_shape *= i;
    }
    return total_shape;
  }
  ~Tensor() {
    if (this->_has_init) {
      free(this->ptr);
    }
  }
  T *ptr;
  Shape shape;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Tensor<T> &obj) {
  os << "[";
  os << obj.ptr[0];
  os << "]";
  return os;
}

typedef std::map<std::string, Tensor<double> *> MapTensor;

#endif
