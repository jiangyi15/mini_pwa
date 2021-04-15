#include <cmath>
#include <iostream>

double fact(int n) {
  double ret = 1.0;
  int m = n / 2;
  for (int i = 1; i <= m; i++) {
    ret *= i;
  }
  return ret;
}

double cg_coeff2(int j1, int j2, int m1, int m2, int j, int m) {
  if (m1 + m2 != m)
    return 0.;
  double a = (j + 1) * fact((j1 + j2 - j)) * fact((j + j1 - j2)) *
             fact((j + j2 - j1)) / fact((j + j1 + j2 + 2));
  a *= fact((j1 + m1)) * fact((j1 - m1)) * fact((j2 + m2)) * fact((j2 - m2)) *
       fact((j + m)) * fact((j - m));
  double b = 0.0;
  int imin = std::max(0, -std::min(j - j2 + m1, j - j1 - m2));
  int imax = std::min(std::min(j1 + j2 - j, j2 + m2), (j1 - m1));
  for (int s = imin; s <= imax; s += 2) {
    double tmp = (s % 2 == 0) ? 1 : -1;
    tmp /= (fact(s)) *
           (fact(j1 + j2 - j - s) * fact(j1 - m1 - s) * fact(j1 + m2 - s) *
            fact(j - j1 - m2 + s) * fact(j - j2 + m1 + s));
    b += tmp;
  }
  return sqrt(a) * b;
}

double cg_coeff(double j1, double j2, double m1, double m2, double j,
                double m) {
  return cg_coeff2(int(j1 * 2), int(j2 * 2), int(m1 * 2), int(m2 * 2),
                   int(j * 2), int(m * 2));
}
