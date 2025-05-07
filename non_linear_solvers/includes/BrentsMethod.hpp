#ifndef BRENTS_METHOD_HPP
#define BRENTS_METHOD_HPP
#include "Function.hpp"

#include <cstddef>

// double secanteCalculation(const Function &, double, double);

// double inverseQuadraticInterpolation(const Function &, double, double,
// double);

bool useBisection(double, double, double, double, double, double, bool);

double BrentsMethod(const FunctionC0 &, double, double, double, size_t);

#endif // BRENTS_METHOD_HPP