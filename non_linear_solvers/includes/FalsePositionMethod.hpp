#ifndef FALSE_POSITION_METHOD_HPP
#define FALSE_POSITION_METHOD_HPP

#include "Function.hpp"
#include <cstddef>

double FalsePositionMethod(const Function &, double, double, double = 0.0001,
                           size_t = 100);

#endif // FALSE_POSITION_METHOD_HPP