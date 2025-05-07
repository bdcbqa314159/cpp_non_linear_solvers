#include "common.hpp"
#include <cassert>
#include <cmath>

double bisection(double a, double b) { return a + (b - a) * 0.5; }

double secante(double a, double b, double f_a, double f_b, double epsilon) {

  double diff = f_b - f_a;
  assert(std::abs(diff) > epsilon);
  double inv_diff_k = 1. / diff;

  return a - f_a * (b - a) / (f_b - f_a);
}