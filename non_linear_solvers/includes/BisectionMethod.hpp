#ifndef BISECTION_METHOD_HPP
#define BISECTION_METHOD_HPP

#include "Function.hpp"
#include <cstddef>

// double bisection(double, double);

double BisectionMethod(const FunctionC0 &, double, double, double = 0.0001,
                       size_t = 100);

#endif // BISECTION_METHOD_HPP