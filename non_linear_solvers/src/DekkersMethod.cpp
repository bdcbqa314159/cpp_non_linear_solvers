#include "DekkersMethod.hpp"
#include <cassert>
#include <algorithm>

double DekkersMethod(const Function& myFunction, double a, double b, double epsilon, size_t maxIterations)
{
    size_t k = 1;
    double a_k = a;
    double b_0 = a, b_1 = b;

    double product_1 = myFunction.Value(a_k)*myFunction.Value(b_1);

    assert(product_1 < 0.);
    if (myFunction.Value(a_k) < myFunction.Value(b_1)) std::swap(a_k, b_1);

    while( (std::abs(myFunction.Value(b_1))>epsilon) && (k<=maxIterations)){

        double diff = myFunction.Value(b_1) - myFunction.Value(b_0);
        double inv_diff = 1./diff;

        double s_k = b_1 - myFunction.Value(b_1)*(b_1 - b_0)*inv_diff;
        double m = a_k + 0.5*(b_1 - a_k);

        b_0 = b_1;

        bool secant_condition = (b_1 > m && m < s_k && s_k < b_1) || (b_1 < m && b_1 < s_k && s_k < m);

        if (secant_condition){
            b_1 = s_k;
        }
        else{
            b_1 = m;
        }

        product_1 = myFunction.Value(a_k)*myFunction.Value(b_1);
        double product_2 = myFunction.Value(b_0)*myFunction.Value(b_1);

        if (product_1 >0. && product_2 < 0.) {a_k = b_0;}

        product_1 = myFunction.Value(a_k)*myFunction.Value(b_1);
        assert(product_1 < 0.);
        if (myFunction.Value(a_k) < myFunction.Value(b_1)) std::swap(a_k, b_1);

    }

    return b_1;
}