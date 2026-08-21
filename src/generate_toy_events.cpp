
#include <generate_toy_events.hpp>

#include <gauss_integrate.hpp>
//ROOT headers
#include <TString.h>
//stdlib headers
#include <stdexcept> 

namespace peak_search
{

void generate_toy_events(Histo1D& hist, Fcn1D* PDF, double N, TRandom3& generator)
{
    //check to make sure N > 0 
    if (N <= 0) {
        throw std::invalid_argument(Form(
            "in <%s>: N <= 0, not allowed! N: %e", 
            __func__, N
        ));
        return; 
    }

    //check to make sure the function passed is not null
    if (!PDF) {
        throw std::invalid_argument(Form(
            "in <%s>: Fcn1D ptr passed (2nd arg) is null!", 
            __func__
        ));
        return; 
    }

    //check to make sure we have at least a few bins
    if (hist.bins.empty()) {
        throw std::invalid_argument(Form(
            "in <%s>: Histo1D passed (1st arg) has no bins!", 
            __func__
        ));
        return; 
    }

    //now, we can do the work we've set out to do. 
    for (auto& bin : hist.bins) {

        double expectation_val = peak_search::gauss_integrate(*PDF, bin.xmin, bin.xmax);

        //check to make sure the expectation valuse isn't negative (can't have a poisson random variable with neg. expecation val!)
        if (expectation_val < 0.) {
            throw std::logic_error(Form(
                "in <%s>: integral of PDF over interval x=[%.6e, x=%.6e] resulted in negative expectation value: %.6e", 
                __func__, bin.xmin, bin.xmax, expectation_val
            ));
            return; 
        }

        bin.N = generator.PoissonD(expectation_val * N); 
    }
    return; 
}

};