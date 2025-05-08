#include "DekkersMethod.hpp"
#include <algorithm>
#include <cassert>

double DekkersMethod(const FunctionC0 &myFunction, double a, double b,
                     double epsilon, size_t maxIterations) {
  size_t k = 1;
  double a_k = a;
  double b_k = a, b_k_plus_1 = b;

  double f_a_k = myFunction(a_k);
  double f_b_k = myFunction(b_k);
  double f_b_k_plus_1 = myFunction(b_k_plus_1);

  double product_1 = f_a_k * f_b_k_plus_1;
  double product_2 = 0.;

  assert(product_1 < 0.);
  if (myFunction(a_k) < myFunction(b_k_plus_1))
    std::swap(a_k, b_k_plus_1);

  f_b_k_plus_1 = myFunction(b_k_plus_1);
  f_a_k = myFunction(a_k);

  while ((std::abs(myFunction(b_k_plus_1)) > epsilon) && (k <= maxIterations)) {

    double s_k = secante(b_k, b_k_plus_1, f_b_k, f_b_k_plus_1, epsilon);
    double m = bisection(a_k, b_k_plus_1);

    b_k = b_k_plus_1;
    f_b_k = myFunction(b_k);

    bool secant_condition = (b_k_plus_1 > m && m < s_k && s_k < b_k_plus_1) ||
                            (b_k_plus_1 < m && b_k_plus_1 < s_k && s_k < m);

    // if no secant condition, use bisection
    if (secant_condition) {
      b_k_plus_1 = s_k;
    } else {
      b_k_plus_1 = m;
    }

    f_b_k_plus_1 = myFunction(b_k_plus_1);

    product_1 = f_a_k * f_b_k_plus_1;
    product_2 = f_b_k * f_b_k_plus_1;

    if (product_1 > 0. && product_2 < 0.) {
      a_k = b_k;
      f_a_k = myFunction(a_k);
    }

    product_1 = f_a_k * f_b_k_plus_1;
    assert(product_1 < 0.);

    if (f_a_k < f_b_k_plus_1) {
      std::swap(a_k, b_k_plus_1);
      f_b_k_plus_1 = myFunction(b_k_plus_1);
      f_a_k = myFunction(a_k);
    }
  }

  return b_k_plus_1;
}