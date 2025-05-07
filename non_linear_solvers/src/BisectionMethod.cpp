#include "BisectionMethod.hpp"
#include "Function.hpp"
#include "common.hpp"
#include <cassert>
#include <cmath>

// double bisection(double a, double b) { return a + (b - a) * 0.5; }

double BisectionMethod(const FunctionC0 &myFunction, double a, double b,
                       double epsilon, size_t maxIterations) {

  assert(a < b);
  size_t k = 1;

  double a_k = a, b_k = b;
  double r_k = a_k;
  double product = 0.;

  double f_ak = myFunction(a_k);
  double f_rk = myFunction(r_k);

  while ((std::abs(f_rk) > epsilon) && (k <= maxIterations)) {

    r_k = bisection(a_k, b_k);
    f_rk = myFunction(r_k);

    product = f_ak * f_rk;

    if (product > 0.) {
      a_k = r_k;
      f_ak = myFunction(a_k);

    } else {
      b_k = r_k;
    }

    k++;
  }

  return r_k;
}