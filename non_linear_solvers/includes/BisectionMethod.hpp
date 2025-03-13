#pragma once
#include <cstddef>
#include "Function.hpp"

double BisectionMethod(const Function &, double, double, double = 0.0001, size_t = 100);