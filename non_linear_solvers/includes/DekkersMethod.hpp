#ifndef DEKKERS_METHOD_HPP
#define DEKKERS_METHOD_HPP

#include "Function.hpp"
#include <cstddef>

double DekkersMethod(const FunctionC0 &, double, double, double = 0.00001,
                     size_t = 100);

#endif // DEKKERS_METHOD_HPP
