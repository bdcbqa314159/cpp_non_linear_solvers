#include "common.hpp"
#include <cassert>
#include <cmath>

double bisection(double a, double b) { return a + (b - a) * 0.5; }

double secante(double a, double b, double f_a, double f_b, double epsilon) {

  double diff = f_b - f_a;
  double inv_diff_k = 1. / diff;

  return a - f_a * (b - a) / (f_b - f_a);
}

double inverseQuadraticInterpolation(const FunctionC0 &myFunction, double a,
                                     double b_0, double b_1) {

  double fa = myFunction(a);
  double fb_0 = myFunction(b_0);
  double fb_1 = myFunction(b_1);

  double P1 = a * fb_0 * fb_1 / ((fa - fb_0) * (fa - fb_1));
  double P2 = b_1 * fa * fb_0 / ((fb_1 - fa) * (fb_1 - fb_0));
  double P3 = b_0 * fa * fb_1 / ((fb_0 - fa) * (fb_0 - fb_1));

  return P1 + P2 + P3;
}

bool useBisection(double a, double b_0, double b_1, double b_2, double s,
                  double epsilon, bool bisection) {

  double delta = epsilon + std::numeric_limits<double>::min();

  bool P1 = bisection && (std::abs(s - b_2) >= 0.5 * std::abs(b_2 - b_1));
  bool P2 = !bisection && (std::abs(s - b_2) >= 0.5 * std::abs(b_1 - b_0));
  bool P3 = bisection && (std::abs(b_2 - b_1) < delta);
  bool P4 = !bisection && (std::abs(b_1 - b_0) < delta);

  return P1 || P2 || P3 || P4;
}
