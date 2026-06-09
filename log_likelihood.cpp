
#include "log_likelihood.hpp"
#include "gauss_integrate.hpp"

//ROOT headers
#include <TAxis.h>

//stdlib headers
#include <stdexcept>
#include <limits> 

namespace peak_search
{

//_______________________________________________________________________________________________
double log_likelihood(TH1D* hist, double (*fcn)(double))
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

        log_likelihood += -expect + ni*std::log(expect) - log_factorial(ni);  
    }   
    return log_likelihood; 
}

//_______________________________________________________________________________________________
double log_likelihood(TH2D* hist, double (*fcn)(double,double))
{   
    if (!hist) {
        throw std::logic_error("in <log_likelihood>: hist ptr is null"); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }

    auto xax = hist->GetXaxis(); 
    auto yax = hist->GetYaxis(); 

    double dx = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins());
    double dy = (yax->GetXmax() - yax->GetXmin())/((double)yax->GetNbins());    

    double log_likelihood = 0.; 

    double x = xax->GetXmin(); 
    for (int bx=1; bx<=xax->GetNbins(); bx++) {
        double y = yax->GetXmin();

        for (int by=1; by<=yax->GetNbins(); by++) {
            double y = yax->GetBinCenter(by);        
        
            double expect = gauss_integrate(fcn, x,x+dx, y,y+dy);

            //number of events in this bin 
            double ni = hist->GetBinContent(bx,by);

            log_likelihood += -expect + ni*std::log(expect) - log_factorial(ni);  
        
            y += dy; 
        }
        x += dx; 
    }   
    return log_likelihood; 
}

//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
};