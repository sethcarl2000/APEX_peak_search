#include "chisquare.hpp"
#include "numbers.hpp"
#include "gauss_integrate.hpp"

#include <Math/ProbFunc.h>
#include <TString.h> 

#include <stdexcept> 

namespace peak_search
{

double chisquare(TH1D* hist, const std::function<double(double)>& fcn)
{
    if (!hist) {
        throw std::logic_error("in <log_likelihood>: hist ptr is null"); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }

    auto xax = hist->GetXaxis();
    const double dx = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()); 

    int DoF = xax->GetNbins(); 
    double chi2 = 0.;

    double x = xax->GetXmin(); 
    for (int bx=1; bx<xax->GetNbins(); bx++) {

        double expect = gauss_integrate(fcn, x,x+dx);
        x += dx; 

        if (expect < 1e-6) {
            throw std::logic_error(Form("in <%s>: encountered bin with non-positive expectation value; this is illegal for chi^2 test."
                " x=%f, f(x)=%f", __func__, x, expect
            ));
            return numbers::nan; 
        }

        //number of events in this bin 
        double ni = hist->GetBinContent(bx);

        double chi = (ni - expect);

        chi2 += chi*chi / expect;
    }   

    return chi2; 
}

};