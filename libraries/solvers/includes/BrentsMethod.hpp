#ifndef BRENTS_METHOD_HPP
#define BRENTS_METHOD_HPP
#include "utils_lib" // IWYU pragma: keep

#include <cstddef>

bool useBisection(double, double, double, double, double, double, bool);

double BrentsMethod(const FunctionC0 &, double, double, double, size_t);

#endif // BRENTS_METHOD_HPP