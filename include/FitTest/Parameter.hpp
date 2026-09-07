#ifndef peak_search_FitTest_Parameter_hpp
#define peak_search_FitTest_Parameter_hpp

#include <numbers.hpp>
// stdlib
#include <stdexcept> 

namespace peak_search
{
namespace FitTest
{

class Parameter {
    unsigned long n_steps;
    double minval, maxval, dx;
public: 
    Parameter(unsigned long n, double min, double max)
        : n_steps{n}, minval{min}, maxval{max}, dx{(max-min)/((double)(n-1))} 
    {
        if (n_steps < 2)
            throw std::invalid_argument("in <FitTest::Parameter(ctor)>: Tried to construct parameter with less than 2 steps, this is not allowed."); 
    }; 

    double get_step(unsigned long i) const {
        if (i >= n_steps) { return numbers::nan; }
        return minval + ((double)i)*dx; 
    }

    unsigned long get_n_steps() const { return n_steps; }
};

}
}



#endif