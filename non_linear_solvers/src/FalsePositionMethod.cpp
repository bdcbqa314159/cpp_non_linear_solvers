#include "FalsePositionMethod.hpp"
#include <cassert>
#include <cmath>

double FalsePositionMethod(const Function& myFunction, double a, double b, double epsilon, size_t maxIterations){

    assert(a<b);

    size_t k = 1;
    double a_k = a, b_k = b;
    double r_k = a_k;

    double diff_k = myFunction.Value(b_k) - myFunction.Value(a_k);

    assert( std::abs(diff_k) > epsilon);
    double inv_diff_k = 1./diff_k;
    

    while( (std::abs(myFunction.Value(r_k))>epsilon) && (k<=maxIterations)){

        r_k = a_k - myFunction.Value(a_k)*(b_k - a_k)*inv_diff_k;
        double product = myFunction.Value(a_k)*myFunction.Value(r_k);

        if (product > 0.){
            a_k = r_k;        
        }
        else{
            b_k = r_k;
        }

        diff_k = myFunction.Value(b_k) - myFunction.Value(a_k);
        assert( std::abs(diff_k) > epsilon);
        inv_diff_k = 1./diff_k;
        k++;
    }

    return r_k;

}