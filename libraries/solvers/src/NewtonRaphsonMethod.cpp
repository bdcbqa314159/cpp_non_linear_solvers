#include "NewtonRaphsonMethod.hpp"
#include <cassert>

double NewtonRaphsonMethod(const FunctionC1 &myFunction, double target,
                           double guess, double epsilon, size_t maxIterations) {

  size_t k = 1;

  double x_prev = guess;
  double f_x_prev = myFunction(x_prev);
  double df_x_prev = myFunction.derivative(x_prev);

  assert(df_x_prev != 0.0);
  double x_current = x_prev - (f_x_prev - target) / df_x_prev;

  while ((std::abs(x_current - x_prev) > epsilon) && (k <= maxIterations)) {

    x_prev = x_current;
    f_x_prev = myFunction(x_prev);
    df_x_prev = myFunction.derivative(x_prev);

    x_current = x_prev - (f_x_prev - target) / df_x_prev;
    k++;
  }

  return x_current;
}