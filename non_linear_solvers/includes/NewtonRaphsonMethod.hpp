#pragma once
#ifndef NEWTONRAPHSON_HPP
#define NEWTONRAPHSON_HPP

#include "Function.hpp"
#include <cstddef>

double NewtonRaphsonMethod(const Function &, double, double, double = 0.0001,
                           size_t = 100);

#endif // NEWTONRAPHSON_HPP