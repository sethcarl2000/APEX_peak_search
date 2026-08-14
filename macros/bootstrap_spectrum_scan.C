
#include "ToyEventGenerator.h"

#include <bininfo.hpp>
#include <copy_subhist.hpp>
#include <fit_exponential_poly.hpp>
#include <fit_exponential_legendre.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <newton_optimizer.hpp> 
#include <fit_parameter.hpp> 
#include <compute_Q0.hpp>
#include <numbers.hpp>
#include <read_model_from_file.hpp>
#include <generate_toy_events.hpp>
//ROOT headers
#include <TRandom3.h> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TAxis.h> 
#include <TGraph.h> 
#include <TLine.h> 
#include <TFile.h>
#include <TList.h> 
#include <TObject.h> 
#include <TPad.h> 
#include <TStyle.h> 
//stdlib headers
#include <cstdio> 
#include <functional> 
#include <iostream> 
#include <cmath> 
#include <thread> 
#include <stdexcept> 
#include <thread> 
#include <mutex> 
#include <sstream> 

#define VERBOSE 1

//#define MAX_THREADS 1

namespace {

    //min & max m-values (taken from accidental spectrum)
    constexpr double m_min = 140, m_max = 280;

    //model for the exponential background
    const std::string path_accidental_model = "data/models/exp_poly_19.dat";
}; 

//this manages our temporary, drawn objects so that each new frame in a gif can have fresh objects. 
class DrawnObjectMgr {
private: 
    TList fDrawnObjects; 
public: 
    
    inline void Draw(TObject* obj, const char* option="") {
        obj->SetBit(kMustCleanup); 
        obj->ResetBit(kCanDelete);
        fDrawnObjects.Add(obj);
        obj->Draw(option);
    }   

    inline void Reset() { fDrawnObjects.Delete(); }
}; 

//get the mass resolution as a fcn of sigma
double GetSigma(double m); 

/// @brief Draw progress bar 
/// @param progress progress amount (in range [0-1])
/// @param n_ticks number of 'ticks' to drawin the progress bar
/// @return formatted progress bar string
std::string progress_bar(double progress, int n_steps=100); 

/// @brief Create a local, deletable TObject
template<typename T> T* local_object_copy(TObject* obj, size_t t)
{
    const char* new_obj_name = Form("%s_t%zi",obj->GetName(),t); 
    auto cpy = dynamic_cast<T*>(obj->Clone(new_obj_name)); 
    cpy->SetBit(kMustCleanup); 
    cpy->ResetBit(kCanDelete);
    cpy->SetDirectory(0); 
    return cpy; 
}


void bootstrap_spectrum_scan(std::string file_path, size_t n_scans=100, std::string histogram_name="h_m")
{   

    //mass ranges to fit
    double min_fit_mass{150}, max_fit_mass{270};
    
    double sigma=1.; 
    int n_steps=400; 

    
    double window_size = sigma*7.; 

    //cx 03.22222222222222222222 
    // - Muon's comment (23 Jul 26)
    std::cout << 
        "getting events..." << std::flush;
    

    peak_search::ExponentialPoly accidental_bg_model({}, m_min, m_max);

    try {

        peak_search::read_model_from_file(path_accidental_model, &accidental_bg_model);

    } catch (const std::exception& e) {

        Error(__func__, "Something went wrong trying to load model from file\n what(): %s", e.what()); 
        return; 
    }
    
    //bin size
    const double dm = 0.5; 
    const size_t n_bins = (m_max - m_min)/dm; 

    peak_search::histo_1D_t total_hist_1d{.bins={}, .xmin=m_min, .xmax=m_max};
    
    total_hist_1d.bins.reserve(n_bins); 

    double m = m_min + dm/2.;
    for (size_t i=0; i<n_bins; i++) {
        
        total_hist_1d.bins.push_back({ m, 0. });
        m += dm; 
    } 

    auto hist = new TH1D("h_test", "Test of toy-event generator", n_bins, m_min, m_max); 

    double n_events_total = 246e6; 

    TRandom3 generator; 
    peak_search::generate_toy_events(total_hist_1d, &accidental_bg_model, n_events_total, generator);

    for (const auto& bin : total_hist_1d.bins) {
        
        hist->Fill( bin.x, bin.N );
    }
    hist->SetBinErrorOption( TH1::kPoisson );

    std::printf("total integral of test-hist: %.0f", hist->Integral());

    auto c_fit = new TCanvas; 
    gStyle->SetOptStat(0); 
    c_fit->DivideRatios(1,2, {1.}, {0.25, 0.75}, 0.00,0.00); 
    
    const double max_signal_events = 40e3; 
    auto hist_S = new TH2D("h_signal", "Best-fit signal parameter '#mu' vs m;signal mass hypothesis (MeV);best-fit #mu", 
        n_steps/4, min_fit_mass, max_fit_mass, 
        100, -max_signal_events, max_signal_events
    ); 

    const double max_significance = 6.;
    auto hist_Z = new TH2D("h_Z", "Significance Z ~ #sqrt{Q0} vs m;signal mass hypothesis (MeV);Significance Z (n. #sigma)", 
        n_steps/4, min_fit_mass, max_fit_mass, 
        100, -max_significance, +max_significance
    ); 

    auto hist_pZ = new TH2D("h_pZ", "p(Q0) vs m;signal mass hypothesis (MeV);p(Q0)", 
        n_steps/8, min_fit_mass, max_fit_mass,
        50, 0., 1.
    ); 

    auto c_fit_hist = c_fit->cd(2);
    c_fit_hist->SetTopMargin(0.);  
    hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum()*1.1); 
    hist->Draw("HIST"); 

    //this histogram will track the dist. of p(Q0). 
    auto hist_pQ0 = new TH1D("h_pQ0", "Dist. of p(Q0);p(Q0);", 50, 0., 1.); 
    auto hist_sqrtQ0 = new TH1D("h_sqrtQ0", "Significance Z = #sqrt{Q0}", 100, -6, 6); 

    std::cout << "done.\n" << std::flush; 

    double xmin{m_min}, xmax{m_max};

    const double dx = (xmax - xmin)/((double)hist->GetXaxis()->GetNbins()); 
    
    //now, get ready to launch threads 
#ifdef MAX_THREADS
    const size_t n_threads = std::min<size_t>( std::thread::hardware_concurrency(), MAX_THREADS ); 
#else 
    const size_t n_threads = std::thread::hardware_concurrency(); 
#endif

    size_t scans_per_thread = n_scans / n_threads; 

    std::vector<std::thread> threads{}; threads.reserve(n_threads); 

    std::mutex read_mutex, write_mutex; 

    size_t scans_scheduled{0}; 
    
    for (size_t t=0; t<n_threads; t++) {

        threads.emplace_back([&]{

            auto t_this = t; 

            //make a copy of the histogram.
            //it may seem dumb to need a mutex for a read operation, but we can't gurantee that the input histogram
            // is 'const' for this operation, thus this headache. 
            read_mutex.lock(); 
            auto data = peak_search::copy_1D_hist(hist); 
            auto h_pQ0_t    = local_object_copy<TH1D>(hist_pQ0, t_this); 
            auto h_sqrtQ0_t = local_object_copy<TH1D>(hist_sqrtQ0, t_this); 
            auto h_t        = local_object_copy<TH1D>(hist, t_this); 
            auto h_S_t      = local_object_copy<TH2D>(hist_S, t_this); 
            auto h_Z_t      = local_object_copy<TH2D>(hist_Z, t_this); 
            auto h_pZ_t     = local_object_copy<TH2D>(hist_pZ, t_this); 
            read_mutex.unlock(); 

            TRandom3 rand_t; 

            //randomly re-sample each datapoint for each scan. 
            
            //scans done by this particular thread
            size_t scans_completed_t {0};

            while (1) {

                read_mutex.lock();
                if (++scans_scheduled > n_scans) { 
                    read_mutex.unlock(); 
                    break; 
                }
                if (VERBOSE == 1) {
                    //std::printf("thread %2zi/%zi performing scan %zi (%zi) (%4.1f%%)\n", t_this+1,n_threads, scans_completed_t+1, scans_scheduled, 100.*((double)scans_scheduled)/((double)n_scans));
                    std::cout << "\r" << progress_bar(((double)scans_scheduled)/((double)n_scans), 100) << std::flush; 
                }
                read_mutex.unlock(); 

                auto h_scan = local_object_copy<TH1D>(h_t, t_this); 

                //________________________________________________________________________________________________________________________________
                auto make_toy_hist = [&rand_t, n_events_total, dm, &accidental_bg_model](double m_lo, double m_hi) {
                    //make sure the range is an even multiple of the number of bins
                    
                    //make sure neither limit is out-of-range
                    if (m_lo < m_min) m_lo = m_min; 
                    if (m_hi > m_max) m_hi = m_max; 

                    //make sure the span of this hist represents an integer number of bins
                    size_t n_bins = (m_hi - m_lo)/dm;   

                    double m_center = (m_hi + m_lo)/2.; 
                    double m_span   = ((double)n_bins)*dm; 

                    m_lo = m_center - m_span/2.; 
                    m_hi = m_center + m_span/2.; 

                    peak_search::histo_1D_t my_hist{ .bins={}, .xmin=m_lo, .xmax=m_hi };

                    my_hist.bins.reserve(n_bins);
                    double m=m_lo + dm/2.;
                    for (int i=0; i<n_bins; i++) {
                        my_hist.bins.push_back({.x=m}); 
                        m += dm; 
                    }

                    //generate the random events
                    peak_search::generate_toy_events(my_hist, &accidental_bg_model, n_events_total, rand_t); 

                    return my_hist; 
                };  
                //_________________________________________________________________________________________________________________________________

                double dm_step = (max_fit_mass - min_fit_mass)/((double)n_steps-1);

                double m = min_fit_mass;

                auto gaussian_fcn   = peak_search::Gauss(0, m, sigma); 
                
                for (int i_step=0; i_step < n_steps; i_step++) {
    
                    gaussian_fcn.Set_x0(m);
                    gaussian_fcn.Set_mu(0.);
                
                    //generate data
                    auto sub_data = make_toy_hist(m - window_size, m + window_size);

                    if (VERBOSE >= 2) {
                        std::cout << "fit m=" << m << "\n";
                        std::cout << "sub histogram xmax=" << m-window_size << ", " << m+window_size << "\n"; 
                        std::cout << "bins:\n"; 
                        for (auto& bin : sub_data.bins) {
                            std::cout << "  " << bin.x << " | " << bin.N << "\n"; 
                        }
                    }

                    auto exp_poly_result = peak_search::fit_exponential_legendre(sub_data, 4);

                    //auto r2 = peak_search::fit_exponential_legendre(sub_data, 4); 

                    if (!exp_poly_result) {
                        Warning(__func__, "m=%.3f fit failed: Poly. background fit failed.", m); 
                        continue; 
                    }

                    auto exp_poly = exp_poly_result.data; 

                    auto fcn_s_plus_b = peak_search::FcnSum(&gaussian_fcn, &exp_poly); 

                    double Q0 = peak_search::compute_Q0(sub_data, fcn_s_plus_b); 

                    if (peak_search::numbers::is_nan(Q0)) {
                        Warning(__func__, "m=%.3f fit failed: computed Q0 is nan", m); 
                        continue; 
                    }

                    double pQ0 = peak_search::compute_Q0_p(Q0);

                    //fill local histogram copies

                    double Z = (Q0>0.?+1.:-1.) * std::sqrt( std::fabs(Q0) ); 

                    // get the 'mu' signal parameter
                    double mu = fcn_s_plus_b.GetParams()[0]; 

                    h_pQ0_t->Fill( pQ0 ); 
                    h_sqrtQ0_t->Fill( Z );

                    h_S_t->Fill( m, mu ); 
                    h_Z_t->Fill( m, Z ); 
                    h_pZ_t->Fill( m, pQ0 ); 

                    m += (max_fit_mass - min_fit_mass)/((double)n_steps-1); 

                } // loop over all choices of 'm' 
                ++scans_completed_t;
            } // loop over all scans
            
            /// Copy contents of one TH1D into another
            auto cpy_hist_1d = [](TH1D* source, TH1D* target) {
                auto xax = source->GetXaxis(); 
                for (int bin=1; bin<=xax->GetNbins(); bin++) {
                    target->Fill( xax->GetBinCenter(bin), source->GetBinContent(bin) );
                }
            };

            /// Copy contents of one TH2D into another 
            auto cpy_hist_2d = [](TH2D* source, TH2D* target) {
                auto xax = source->GetXaxis(); 
                auto yax = source->GetYaxis(); 
                for (int bx=1; bx<=xax->GetNbins(); bx++) {
                    for (int by=1; by<=yax->GetNbins(); by++) {
                        target->Fill( xax->GetBinCenter(bx), yax->GetBinCenter(by), source->GetBinContent(bx, by) );
                    }
                }    
            };

            write_mutex.lock();
            //if (VERBOSE >= 1) std::printf("thread %2zi/%zi done with all scans.\n", t+1, n_threads); 
            cpy_hist_1d(h_pQ0_t, hist_pQ0);
            cpy_hist_1d(h_sqrtQ0_t, hist_sqrtQ0); 
            cpy_hist_2d(h_S_t, hist_S);
            cpy_hist_2d(h_Z_t, hist_Z); 
            cpy_hist_2d(h_pZ_t, hist_pZ);
            write_mutex.unlock(); 
            
            delete h_pQ0_t; 
            delete h_sqrtQ0_t;
            delete h_t; 
            delete h_S_t; 
            delete h_Z_t; 
        }); 
    }
    for (auto& t : threads) t.join(); 

    if (VERBOSE==1) std::cout << "\n"; 

    auto reset_bin_errors = [](TH1D* h) {
        auto xax = h->GetXaxis(); 
        for (int bin=1; bin<=xax->GetNbins(); bin++) { h->SetBinError(bin, std::sqrt(h->GetBinContent(bin))); }
    };

    new TCanvas; 
    reset_bin_errors(hist_pQ0); 
    hist_pQ0->SetMaximum( hist_pQ0->GetMaximum()*1.2 ); 
    hist_pQ0->SetMinimum( 0. ); 
    hist_pQ0->Draw("HIST, E");

    new TCanvas; 
    reset_bin_errors(hist_sqrtQ0); 
    hist_sqrtQ0->Draw("HIST, E"); 

    new TCanvas; 
    hist_S->Draw("col"); 

    new TCanvas; 
    hist_Z->Draw("col");

    new TCanvas; 
    hist_pZ->Draw("col");

    return; 
}


double GetSigma(double m) {
    return 1.; 
}

std::string progress_bar(double progress, int n_ticks) {
    std::ostringstream oss; 
    oss  << "["; 
    double n_steps_d = (double)n_ticks;
    int n_ticks_full = ((double)n_ticks)*progress;

    for (int i=0; i<n_ticks; i++) { oss << (i<=n_ticks_full ? "=" : " "); }

    oss << Form("]    %4.1f%%", progress*100.); 
    return oss.str(); 
}
