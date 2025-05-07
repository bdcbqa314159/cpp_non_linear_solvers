#include "BrentsMethod.hpp"
#include "common.hpp"
#include <cassert>
#include <cmath>

double secanteCalculation(const Function &myFunction, double a, double b) {
  return b - myFunction.Value(b) * (b - a) /
                 (myFunction.Value(b) - myFunction.Value(a));
}

double inverseQuadraticInterpolation(const Function &myFunction, double a,
                                     double b_0, double b_1) {

  double fa = myFunction.Value(a);
  double fb_0 = myFunction.Value(b_0);
  double fb_1 = myFunction.Value(b_1);

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

double BrentsMethod(const Function &myFunction, double a, double b,
                    double epsilon, size_t maxIterations) {

  double product = myFunction.Value(a) * myFunction.Value(b);

  assert(product < 0.);
  if (myFunction.Value(a) < myFunction.Value(b))
    std::swap(a, b);

  size_t k = 1;
  double a_k = a;
  double b_0 = a, b_1 = a, b_2 = b;
  double s = std::numeric_limits<double>::max();
  double fs = std::numeric_limits<double>::max();

  bool bisection = true;

  while ((std::abs(myFunction.Value(b_1)) > epsilon) &&
         std::abs(fs) > epsilon && std::abs(b_2 - a_k) > epsilon &&
         (k <= maxIterations)) {

    bool inverse_quadratic_interpolation_condition =
        (myFunction.Value(a_k) != myFunction.Value(b_1) &&
         myFunction.Value(b_2) != myFunction.Value(b_1));

    if (inverse_quadratic_interpolation_condition) {
      s = inverseQuadraticInterpolation(myFunction, a_k, b_1, b_2);
    } else {
      s = secanteCalculation(myFunction, a_k, b_2);
    }

    bool bisection_condition =
        useBisection(a_k, b_0, b_1, b_2, s, epsilon, bisection);

    if (bisection_condition) {
      s = a_k + 0.5 * (b_2 - a_k);
      bisection = true;
    } else {
      bisection = false;
    }

    fs = myFunction.Value(s);
    b_0 = b_1;
    b_1 = b_2;

    product = myFunction.Value(a_k) * fs;
    if (product < 0.)
      b_2 = s;
    else
      a_k = s;

    product = myFunction.Value(a_k) * myFunction.Value(b_2);
    assert(product < 0.);
    if (myFunction.Value(a_k) < myFunction.Value(b_2))
      std::swap(a_k, b_2);
  }

  double fb = myFunction.Value(b_2);

  return fb < fs ? b : s;
}

double newBrentsMethod(const FunctionC0 &myFunction, double a, double b,
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

  while ((std::abs(myFunction.Value(b_1)) > epsilon) &&
         std::abs(fs) > epsilon && std::abs(b_2 - a_k) > epsilon &&
         (k <= maxIterations)) {

    bool inverse_quadratic_interpolation_condition =
        (myFunction.Value(a_k) != myFunction.Value(b_1) &&
         myFunction.Value(b_2) != myFunction.Value(b_1));

    if (inverse_quadratic_interpolation_condition) {
      s = inverseQuadraticInterpolation(myFunction, a_k, b_1, b_2);
    } else {
      s = secanteCalculation(myFunction, a_k, b_2);
    }

    bool bisection_condition =
        useBisection(a_k, b_0, b_1, b_2, s, epsilon, bisection);

    if (bisection_condition) {
      s = a_k + 0.5 * (b_2 - a_k);
      bisection = true;
    } else {
      bisection = false;
    }

    fs = myFunction.Value(s);
    b_0 = b_1;
    b_1 = b_2;

    product = myFunction.Value(a_k) * fs;
    if (product < 0.)
      b_2 = s;
    else
      a_k = s;

    product = myFunction.Value(a_k) * myFunction.Value(b_2);
    assert(product < 0.);
    if (myFunction.Value(a_k) < myFunction.Value(b_2))
      std::swap(a_k, b_2);
  }

  double fb = myFunction.Value(b_2);

  return fb < fs ? b : s;
}