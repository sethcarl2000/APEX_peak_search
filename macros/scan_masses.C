#ifndef scan_masses_C
#define scan_masses_C

#include "../chisquare.hpp"
#include "gauss_integrate.hpp"
#include "compute_Q0.hpp"
#include "fit_exponential_poly.hpp"
#include "ExponentialPoly.hpp"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h> 
#include <TRandom3.h> 
#include <Math/ProbFunc.h>
#include <Math/QuantFuncMathCore.h>
#include <Math/SpecFuncMathCore.h>
#include <Math/PdfFuncMathCore.h>
#include <TF1.h> 
#include <TAxis.h> 
#include <TError.h> 

#include <cmath> 
#include <cstdio> 
#include <functional>
#include <vector>  

//get the width of the signal mass peak
double mass_width(double m) {
    return (( 0.7 - 1.0 )/(250. - 140.)) * (m - 140.) + 1.; 
};

double normal_pdf(double x, double sigma) {
    x = x/sigma; 
    return std::exp( - x*x/2. ) / ( 2.50662827463 * sigma ); 
}

struct m_and_Q0_t { double m, Q0, p_Q0, chi2, p_chi2; };

/// @brief Scan mass hypothesis over a given range, return a list of 
/// @param hist_m_p histogram, where x-axis is mass, and y-axis is Pa (accidental probability)
/// @param n_tests //number of mass-tests to perform 
/// @param m_min //minimum mass to test
/// @param m_max //maximum mass to test 
/// @param max_Pa //cut on Pa (all events with a Pa smaller than this value are kept.)
/// @return vector of mass hypothesis, and associated Q0-values.
std::vector<m_and_Q0_t> scan_masses(TH2D* hist_m_p, int n_tests, double m_min, double m_max, double max_Pa, int exp_bg_order=6)
{
    if (!hist_m_p) return {}; 

    std::vector<m_and_Q0_t> experiments; experiments.reserve(n_tests); 

    auto xax = hist_m_p->GetXaxis();
    auto yax = hist_m_p->GetYaxis(); 
    
    const int n_bins_x = xax->GetNbins(); 

    //first, flatten the TH2D into a 1D histogram 
    TH1D *mass_hist = (TH1D*)hist_m_p->ProjectionX(Form("%s_proj",hist_m_p->GetTitle()), 0, yax->FindBin(max_Pa)+1);
    mass_hist->SetDirectory(0); 

    double window_scale = 6.; 
    
    double dx = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins());
    
    double avg_p_accidental=0.; 

    //__________________________________________________________________________________________________________
    auto create_sub_hist = [dx,xax,mass_hist](double m0, double m1) {
        m0 = std::max( xax->GetXmin(), m0 ); 
        m1 = std::min( xax->GetXmax(), m1 ); 
        
        int min_bin = xax->FindBin(m0);
        int max_bin = xax->FindBin(m1);

        m0 = xax->GetBinCenter(min_bin)-dx/2.;
        m1 = xax->GetBinCenter(max_bin)+dx/2.;
        
        double center = (m1 + m0)/2.;
        double scale  = (m1 - m0)/2.;

        peak_search::histo_1D_t data;
        data.xmin = -1.; 
        data.xmax = +1.; 
        data.bins.reserve(max_bin-min_bin+1);  

        for (int ix=min_bin; ix<=max_bin; ix++) {
            data.bins.emplace_back( (xax->GetBinCenter(ix)-center)/scale, mass_hist->GetBinContent(ix) );
        }
        return data; 
    };
    //__________________________________________________________________________________________________________

    double m = m_min; 
    double dm = (m_max - m_min)/((double)n_tests-1);
    for (int i=0; i<n_tests; i++) {

        //create a sub-histogram
        double window_size = window_scale * mass_width(m); 

        auto data = create_sub_hist( m - window_size, m + window_size ); 

        auto fit_result = peak_search::fit_exponential_poly(data, exp_bg_order);         
        
        if (!fit_result) {
            continue;  
        } 

        auto poly = (peak_search::ExponentialPoly)fit_result; 

        double chi2 = peak_search::chisquare(data, poly); 
        double p_chi2 = peak_search::chisquare_p(chi2, data.bins.size()); 

        auto nuissance_params = poly.coeffs; 

        double sigma_m = mass_width(m); 
        auto Q0_wrapper = (std::function<double(double,double,const double*)>)[m, sigma_m, exp_bg_order](double x, double mu, const double* par){
            return peak_search::ExponentialPoly::Eval(x, par, exp_bg_order) + mu*normal_pdf(x, sigma_m); 
        };

        double Q0 = peak_search::compute_Q0(data, Q0_wrapper, nuissance_params); 
         
        double p_Q0 = peak_search::compute_Q0_p(Q0);
        //std::printf(" pQ0: %f Q0: %f\n", p_Q0, Q0); 
        experiments.emplace_back( m, Q0, p_Q0, chi2, p_chi2 );  
        m += dm; 
    }
    delete mass_hist; 

    return experiments; 
}

#endif 