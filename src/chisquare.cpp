#include "chisquare.hpp"
#include "numbers.hpp"
#include "gauss_integrate.hpp"
#include <make_histogram_copy.hpp>
// ROOT headers
#include <Math/ProbFunc.h>
#include <Math/QuantFuncMathCore.h>
#include <TString.h> 
// stdlib headers
#include <stdexcept> 

namespace peak_search
{

//_______________________________________________________________________________________-
double chisquare(TH1D* hist, const Fcn1D& fcn)
{
    if (!hist) {
        throw std::logic_error("in <chisquare>: hist ptr is null"); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }

    Histo1D data = make_histogram_copy(hist); 

    return chisquare(data, fcn); 
}
//_______________________________________________________________________________________-
double chisquare(const Histo1D& data, const Fcn1D& fcn)
{
    int DoF = data.bins.size(); 
    double chi2 = 0.;
 
    for (const auto& bin : data.bins) {

        double expect = gauss_integrate(fcn, bin.xmin, bin.xmax);

        if (expect < 1e-6) {
            throw std::logic_error(Form("in <%s>: encountered bin with non-positive expectation value; this is illegal for chi^2 test."
                " x=%f, f(x)=%f", __func__, (bin.xmax + bin.xmin)/2., expect
            ));
            return numbers::nan; 
        }

        //number of events in this bin 
        double ni = bin.N;

        double chi = (ni - expect);

        chi2 += chi*chi / expect;
    }   

    return chi2; 
}
//_______________________________________________________________________________________-
double chisquare_p(double chi2, int n_bins)
{
    //return ROOT::Math::gamma_cdf_c(chi2/2., ((double)n_bins)/2., 1.); 
    return ROOT::Math::chisquared_cdf(chi2, n_bins);
}

};