#include "Function.hpp"

FunctionC0::FunctionC0(const std::function<double(double)> &f) : f(f) {}
double FunctionC0::operator()(double x) const { return f(x); }

FunctionC1::FunctionC1(const std::function<double(double)> &_f,
                       const std::function<double(double)> &_df)
    : FunctionC0(_f), df(_df) {}

double FunctionC1::derivative(double x) const { return df(x); }
