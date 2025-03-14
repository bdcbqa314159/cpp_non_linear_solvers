#pragma once
#include <cstddef>
#include "Function.hpp"

double DekkersMethod(const Function& , double , double , double = 0.00001, size_t = 100);
