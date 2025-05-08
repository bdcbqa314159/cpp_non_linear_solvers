#include "BrentsMethod.hpp"
#include <cassert>
#include <cmath>

double BrentsMethod(const FunctionC0 &myFunction, double a, double b,
                    double epsilon, size_t maxIterations) {

  double f_a = myFunction(a);
  double f_b = myFunction(b);
  double product = f_a * f_b;

  assert(product < 0.);

  if (f_a < f_b)
    std::swap(a, b);

  size_t k = 1;
  double a_k = a;
  double b_0 = a, b_1 = a, b_2 = b;
  double s = std::numeric_limits<double>::max();
  double fs = std::numeric_limits<double>::max();

  bool bisection = true;
  double f_b_1 = myFunction(b_1);
  double f_b_2 = myFunction(b_2);
  double f_a_k = myFunction(a_k);

  while ((std::abs(f_b_1) > epsilon) && std::abs(fs) > epsilon &&
         std::abs(b_2 - a_k) > epsilon && (k <= maxIterations)) {

    bool inverse_quadratic_interpolation_condition =
        (f_a_k != f_b_1 && f_b_2 != f_b_1);

    if (inverse_quadratic_interpolation_condition) {
      s = inverseQuadraticInterpolation(myFunction, a_k, b_1, b_2);
    } else {
      s = secante(a_k, b_2, f_a_k, f_b_2, epsilon);
    }

    bool bisection_condition =
        useBisection(a_k, b_0, b_1, b_2, s, epsilon, bisection);

    if (bisection_condition) {
      s = a_k + 0.5 * (b_2 - a_k);
      bisection = true;
    } else {
      bisection = false;
    }

    fs = myFunction(s);
    b_0 = b_1;
    b_1 = b_2;

    f_b_1 = myFunction(b_1);

    product = f_a_k * fs;
    if (product < 0.) {
      b_2 = s;
      f_b_2 = myFunction(b_2);
    } else {
      a_k = s;
      f_a_k = myFunction(a_k);
    }

    product = f_a_k * f_b_2;
    assert(product < 0.);
    if (f_a_k < f_b_2)
      std::swap(a_k, b_2);
  }

  double fb = f_b_2;

  return fb < fs ? b : s;
}