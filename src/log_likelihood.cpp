
#include "log_likelihood.hpp"
#include "gauss_integrate.hpp"
#include "numbers.hpp"

//ROOT headers
#include <TAxis.h>

//stdlib headers
#include <stdexcept>
#include <limits> 

namespace peak_search
{

//_______________________________________________________________________________________________
double negative_log_likelihood(const Histo1D& hist, const std::function<double(double)>& fcn)
{   
#ifdef DEBUG
    std::printf("in <negative_log_likelihood>\n"); 
#endif

    double log_likelihood = 0.; 

    for (const auto& bin : hist.bins) {
        double expect = gauss_integrate(fcn, bin.xmin, bin.xmax);

        //number of events in this bin 
        double ni = bin.N;

        log_likelihood += -expect + ni*std::log(expect) - numbers::log_factorial(ni);  
    }   
    return -log_likelihood; 
}
//_______________________________________________________________________________________________
double modified_nll(const Histo1D& hist, const std::function<double(double)>& fcn)
{   
#ifdef DEBUG
    std::printf("in <negative_log_likelihood>\n"); 
#endif
    
    double NLL = 0.; 

    for (const auto& bin : hist.bins) {
        
        double expect = gauss_integrate(fcn, bin.xmin, bin.xmax);

        //number of events in this bin 
        double ni = bin.N;

        NLL += expect - ni*std::log(expect);
    }
    return NLL; 
}
//_______________________________________________________________________________________________
double log_likelihood(TH1D* hist, const std::function<double(double)>& fcn)
{   
    if (!hist) {
        throw std::logic_error("in <log_likelihood>: hist ptr is null"); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }

    auto xax = hist->GetXaxis();
    const double dx = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()); 

    double log_likelihood = 0.; 

    double x = xax->GetXmin(); 
    for (int bx=1; bx<xax->GetNbins(); bx++) {

        double expect = gauss_integrate(fcn, x,x+dx);
        x += dx; 

        //number of events in this bin 
        double ni = hist->GetBinContent(bx);

        log_likelihood += -expect + ni*std::log(expect) - numbers::log_factorial(ni);  
    }   
    return log_likelihood; 
}
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
};