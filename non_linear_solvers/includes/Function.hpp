#ifndef FUNCTION_HPP
#define FUNCTION_HPP

class Function {

public:
  virtual double Value(double) const = 0;
  virtual double Derivative(double) const = 0;
};

#endif // FUNCTION_HPP