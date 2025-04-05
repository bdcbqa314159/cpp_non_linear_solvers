#include "NewtonRaphsonMethod.hpp"
#include <cassert>

double NewtonRaphsonMethod(const Function &myFunction, double target,
                           double guess, double epsilon, size_t maxIterations) {

  size_t k = 1;

  double x_prev = guess;
  double x_next = x_prev - (myFunction.Value(x_prev) - target) /
                               myFunction.Derivative(x_prev);

  while (abs(x_next - x_prev) > epsilon && (k <= maxIterations)) {
    x_prev = x_next;
    x_next = x_prev - (myFunction.Value(x_prev) - target) /
                          myFunction.Derivative(x_prev);
    k++;
  }

  return x_next;
}