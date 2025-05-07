#pragma once
#ifndef NON_LINEAR_SOLVERS_COMMON_HPP
#define NON_LINEAR_SOLVERS_COMMON_HPP

#include "Function.hpp"

double bisection(double, double);

double secante(double, double, double, double, double);

double inverseQuadraticInterpolation(const FunctionC0 &, double, double,
                                     double);

bool useBisection(double, double, double, double, double, double, bool);

#endif // NON_LINEAR_SOLVERS_COMMON_HPP