#include "chisquare.hpp"
#include "numbers.hpp"
#include "gauss_integrate.hpp"

#include <Math/ProbFunc.h>
#include <Math/QuantFuncMathCore.h>
#include <TString.h> 

#include <stdexcept> 

namespace peak_search
{

//_______________________________________________________________________________________-
double chisquare(TH1D* hist, const std::function<double(double)>& fcn)
{
    if (!hist) {
        throw std::logic_error("in <log_likelihood>: hist ptr is null"); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }

    auto data = copy_1D_hist(hist); 

    return chisquare(data, fcn); 
}
//_______________________________________________________________________________________-
double chisquare(histo_1D_t data, const std::function<double(double)>& fcn)
{
    const double dx = (data.xmax - data.xmin)/((double)data.bins.size()); 

    int DoF = data.bins.size(); 
    double chi2 = 0.;

    double x = data.xmin; 
    for (int bx=0; bx<data.bins.size(); bx++) {

        double expect = gauss_integrate(fcn, x,x+dx);
        x += dx; 

        if (expect < 1e-6) {
            throw std::logic_error(Form("in <%s>: encountered bin with non-positive expectation value; this is illegal for chi^2 test."
                " x=%f, f(x)=%f", __func__, x, expect
            ));
            return numbers::nan; 
        }

        //number of events in this bin 
        double ni = data.bins[bx].N;

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