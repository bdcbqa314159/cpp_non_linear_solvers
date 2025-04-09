#include "FalsePositionMethod.hpp"
#include "Function.hpp"
#include <cassert>
#include <cmath>

double secante(double a, double b, double f_a, double f_b, double epsilon) {

  double diff = f_b - f_a;
  assert(std::abs(diff) > epsilon);
  double inv_diff_k = 1. / diff;

  return a - f_a * (b - a) / (f_b - f_a);
}

double FalsePositionMethod(const FunctionC0 &myFunction, double a, double b,
                           double epsilon, size_t maxIterations) {

  assert(a < b);

  size_t k = 1;
  double a_k = a, b_k = b;
  double r_k = a_k;

  double f_rk = myFunction(r_k);
  double f_ak = myFunction(a_k);
  double f_bk = myFunction(b_k);

  double diff_k = f_bk - f_ak;
  assert(std::abs(diff_k) > epsilon);
  double inv_diff_k = 1. / diff_k;
  double product = 0.;

  while ((std::abs(f_rk) > epsilon) && (k <= maxIterations)) {

    r_k = secante(a_k, b_k, f_ak, f_bk, epsilon);
    f_rk = myFunction(r_k);
    product = f_ak * f_rk;

    if (product > 0.) {
      a_k = r_k;
      f_ak = myFunction(a_k);
    } else {
      b_k = r_k;
      f_bk = myFunction(b_k);
    }

    diff_k = f_bk - f_ak;
    assert(std::abs(diff_k) > epsilon);
    inv_diff_k = 1. / diff_k;
    k++;
  }

  return r_k;
}