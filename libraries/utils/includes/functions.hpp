#pragma once
#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <functional>

class Function {

public:
  virtual double Value(double) const = 0;
  virtual double Derivative(double) const = 0;
};

class FunctionC0 {
protected:
  std::function<double(double)> f;

public:
  FunctionC0(const std::function<double(double)> &);
  double operator()(double x) const;
};

class FunctionC1 : public FunctionC0 {
private:
  std::function<double(double)> df;

public:
  FunctionC1(const std::function<double(double)> &,
             const std::function<double(double)> &);

  double derivative(double x) const;
};

#endif // FUNCTION_HPP