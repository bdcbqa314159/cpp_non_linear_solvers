#include "BisectionMethod.hpp"
#include <cassert>
#include <cmath>

double BisectionMethod(const Function& myFunction, double a, double b, double epsilon, size_t maxIterations){

    assert(a<b);
    size_t k = 1;

    double a_k = a, b_k = b;
    double r_k = a_k;

    while( (std::abs(myFunction.Value(r_k))>epsilon) && (k<=maxIterations)){

        r_k = a_k + (b_k - a_k)*0.5;
        double product = myFunction.Value(a_k)*myFunction.Value(r_k);

        if (product > 0.){
            a_k = r_k;        
        }
        else{
            b_k = r_k;
        }

        k++;
    }

    return r_k;

}