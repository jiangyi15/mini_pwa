#include "dfunction.h"

double factorial2(int x) {
  double a = 1.0;
  while (x > 0) {
    a *= x--;
  }
  return a;
}

Tensor<double>
small_d_weight(int j) { //  the prefactor in the d-function of beta
  /*
  For a certain j, the weight coefficient with index (:math:`m_1,m_2,l`) is
  :math:`w^{(j,m_1,m_2)}_{l} =
  (-1)^{m_1-m_2+k}\\frac{\\sqrt{(j+m_1)!(j-m_1)!(j+m_2)!(j-m_2)!}}{(j-m_1-k)!(j+m_2-k)!(m_1-m_2+k)!k!}`,
  and :math:`l` is an integer ranging from 0 to :math:`2j`.

  :param j: Integer :math:`2j` in the formula???
  :return: Of the shape (**j** +1, **j** +1, **j** +1). The indices correspond
  to (:math:`l,m_1,m_2`)
  **/
  auto ret = Tensor<double>(
      {(unsigned int)j + 1, (unsigned int)j + 1, (unsigned int)j + 1});
  auto size = j + 1;

  auto f = [](int x) { return factorial2(x >> 1); };

  int m, n, k;
  for (m = -j; m < j + 1; m += 2) {
    for (n = -j; n < j + 1; n += 2) {
      for (k = max(0, n - m); k < min(j - m, j + n) + 1; k += 2) {
        auto l = (2 * k + (m - n)) / 2;
        auto sign = pow(-1, (k + m - n) / 2);
        auto tmp = sign * sqrt(1.0 * f(j + m) * f(j - m) * f(j + n) * f(j - n));
        tmp /= f(j - m - k) * f(j + n - k) * f(k + m - n) * f(k);
        ret.ptr[l * size * size + (m + j) / 2 * size + (n + j) / 2] = tmp;
      }
    }
  }
  return ret;
}

Tensor<double> small_d_function(int j, Tensor<double> *weight,
                                Tensor<double> *theta) {
  /*
   */
  auto n = theta->shape[0];
  auto size = (unsigned int)j + 1;
  Tensor<double> sincos({size, n});
  Tensor<double> ret({n, size, size});
  for (int i = 0; i < j + 1; i++) {
    for (int k = 0; k < n; k++) {
      sincos.ptr[i * n + k] =
          pow(sin(theta->ptr[k] / 2), i) * pow(cos(theta->ptr[k] / 2), j - i);
    }
  }

  for (int i = 0; i < size; i++) {
    for (int k = 0; k < size; k++) {
      for (int idx = 0; idx < n; idx++) {
        ret.ptr[idx * size * size + i * size + k] = 0;
        for (int l = 0; l < j + 1; l++) {
          ret.ptr[idx * size * size + i * size + k] +=
              sincos.ptr[l * n + idx] *
              weight->ptr[l * size * size + i * size + k];
        }
      }
    }
  }
  return ret;
}
