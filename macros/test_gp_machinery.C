
#include <GP.hpp>
#include <Histo1D.hpp> 
#include <make_histogram_copy.hpp>
#include <numbers.hpp>
// NLopt (nonlinear optimization lib)
#include <nlopt.hpp> 
// Eigen 
#include <eigen3/Eigen/Core>
// Minuit 
#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
// ROOT
#include <TH1D.h> 
#include <TH2D.h> 
#include <TAxis.h> 
#include <TRandom3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h> 
#include <TError.h> 
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
// stdlib 
#include <thread> 
#include <vector> 
#include <cmath> 
#include <functional> 
#include <mutex> 
#include <memory> 

#define DEBUG

/// @brief The objective function to be minimized by 'NLopt' routines 
/// @param n size of parameter-space (size of 'params' array)
/// @param params array of parameters to evaluate at a given point 
/// @param grad gradient at current point (ignore if gradient cannot be computed)
/// @param data arb. data needed by optimization function 
/// @return loss function evaluated at point specified by 'params' 
using nlopt_objective_fcn = double(*)(unsigned, const double*, double*, void*); 

/// @brief Blind window fit alg. 
/// @param n number of hyper-params for gaussian process in question
/// @param params hyperparameters to evaluate 
/// @param grad gradient (ignored) 
/// @param hist data to fit (Histo1D)
/// @return average chi-square fit to blind window using GP with given hyperparams
double blind_window_fit(unsigned n, const double* params, double* grad, void* hist);


double marginal_nll(unsigned n, const double* params, double* grad, void* data); 


// generate data to be added to histogram 
struct blind_window_fitdata {
    // data
    std::vector<peak_search::GP::Point> data; 
    // Kernel function
    peak_search::GP::Kernel kernel{}; 
    // blind window size (number of bins)
    int blind_window_size;
    // minimum and maximum bins to fit  
    int buffer;
};

double gen_event(double min, double max, TRandom3&);

void test_gp_machinery(std::string path_file="data/hist-accidental.root", std::string hist_name="h_m")
{

    
    using namespace peak_search;


    TH1D* hist; 
    TFile* infile; 

    if (path_file.empty()) {

        //first, generate the toy histogram 
        TRandom3 mrand; 

        double xmin{-5.}, xmax{+5.};
        auto hist = new TH1D("h_data", "Toy data to be fit (un-normalized);x;counts", 50, -5., +5.); 

        unsigned long n_events = 100e3; 

        for (unsigned long i=0; i<n_events; i++) hist->Fill( gen_event(xmin, xmax, mrand) ); 
    } else {

        // open the histogram from the file
        infile = new TFile(path_file.c_str(), "READ"); 

        if (!infile->IsOpen()) {
            Error(__func__, "TFile could not be opened: '%s'", path_file.c_str()); 
            return; 
        }

        hist = infile->Get<TH1D>(hist_name.c_str()); 
        if (!hist) { 
            Error(__func__, "Could not find hist '%s' in TFile '%s'", hist_name.c_str(), path_file.c_str()); 
            return; 
        }
    }

    new TCanvas;
    gPad->SetLogy(1); 
    hist->DrawCopy("E"); 

    //normalize the data: 

    GP::TrainParameters params; 


    params.data = make_histogram_copy(hist); 

    std::vector<double> pars{std::log(7.), std::log(2.)}; 

    params.kernel = GP::Kernel{
        [](double x1, double x2, const std::vector<double>& params){
            
            if (params.size() != 2) return peak_search::numbers::nan; 
            double amplitude{std::exp(params[0])}, length_scale{std::exp(params[1])}; 
            double arg = (x1 - x2)/length_scale; 
            return amplitude * std::exp( -0.5*arg*arg ); 
        }, pars
    };
    
    params.params_bounds = {}; 

    params.draw_fit = true; 

    auto computer = std::make_unique<GP::Computer>(); 

    GP::Train(params, computer.get()); 

    if (!params) {
        Error(__func__, "Something went wrong with the GP training. what(): %s", params.status.c_str()); 
    }

    // scan 2d parameter space 

    auto all_points = GP::Normalize_points(params.data); 
    auto nll_computer = std::make_unique<GP::LikelihoodComputer>(all_points, params.kernel); 

    int n_bins = 200; 

    auto hist_param_space = new TH2D("h_data", "Parameter space;log A;log Length;- log L",  
        n_bins, std::log(0.002), std::log(1e6), 
        n_bins, std::log(0.0045), std::log(2.8)
    ); 

    double max_nll = 2.5e3; 

    auto xax = hist_param_space->GetXaxis();
    auto yax = hist_param_space->GetYaxis(); 

    struct nll_pt { double x, y, val; };
    std::vector<nll_pt> nll_pts; nll_pts.reserve(n_bins*n_bins); 
    for (int ix=1; ix<=n_bins; ix++)
        for (int iy=1; iy<=n_bins; iy++)
            nll_pts.emplace_back( xax->GetBinCenter(ix), yax->GetBinCenter(iy), 0. ); 


    const size_t total_points = nll_pts.size(); 

    const size_t n_threads = std::thread::hardware_concurrency(); 

    size_t tasks_per_thread = total_points/n_threads; 

    std::vector<std::thread> threads; threads.reserve(n_threads); 

    std::printf("total points: %zi\n", total_points); 

    size_t start=0; 
    for (size_t t=0; t<n_threads; t++) {

        size_t end = start + tasks_per_thread + (total_points % n_threads > t ? 1 : 0);

        const auto cptr_const = nll_computer.get(); 

        threads.emplace_back([&nll_pts,t,cptr_const,start,end,max_nll]{

            //make a local copy of the likelihood computer
            GP::LikelihoodComputer my_cptr(*cptr_const); 

            std::vector<double> my_pars(2, 0.); 

            for (size_t ti=start; ti<end; ti++) {
            
                auto& point = nll_pts[ti];         
                my_pars[0] = point.x; 
                my_pars[1] = point.y; 
                double nll = my_cptr.ComputeNLL(my_pars.data());   
                if (peak_search::numbers::is_nan(nll) || nll > max_nll) {
                    point.val = 0.; 
                } else { 
                    point.val = nll; 
                }
            }
        }); 

        std::printf("thread %zi/%zi, range=[%zi, %zi]\n", t, n_threads-1, start, end); 
        start = end; 
    }

    for (auto& thread : threads) thread.join(); 

    auto canv = new TCanvas; 
    canv->SetRightMargin(0.15); 
    for (const auto& pt : nll_pts) hist_param_space->Fill( pt.x, pt.y, pt.val ); 
    hist_param_space->Draw("colz"); 
    return; 
}

//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
double gen_event(double min, double max, TRandom3& gen)
{
    double x; 
    do { x = gen.Gaus()*2. - 1.; x += 0.2*std::sin(x*3.1415926536/5.); } while (x > max || x < min); 
    return x; 
}