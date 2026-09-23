
#include <GP/Point.hpp>
#include <GP/Compute.hpp>
#include <Histo1D.hpp> 
#include <make_histogram_copy.hpp>
#include <numbers.hpp>
// NLopt (nonlinear optimization lib)
#include <nlopt.hpp> 
// ROOT
#include <TH1D.h> 
#include <TRandom3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h> 
#include <TError.h> 
#include <TCanvas.h>
// stdlib 
#include <vector> 
#include <cmath> 
#include <functional> 
#include <mutex> 

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

struct Kernel {
    std::function<double(double,double,const std::vector<double>&)> fcn{
        [](double,double,const std::vector<double>&){ return peak_search::numbers::nan; }
    }; 
    std::vector<double> params{}; 

    double operator()(double x1, double x2) const { return fcn(x1,x2,params); }
}; 
static_assert(std::is_default_constructible_v<Kernel>); 

// generate data to be added to histogram 
struct blind_window_fitdata {
    // data
    std::vector<peak_search::GP::Point> data; 
    // Kernel function
    Kernel kernel{}; 
    // blind window size (number of bins)
    int blind_window_size;
    // minimum and maximum bins to fit  
    int buffer;
};

double gen_event(double min, double max, TRandom3&);

void test_train_GP()
{
    using namespace peak_search;

    //first, generate the toy histogram 
    TRandom3 mrand; 

    double xmin{-5.}, xmax{+5.};
    auto hist = new TH1D("h_data", "Toy data to be fit (un-normalized);x;counts", 50, -5., +5.); 

    unsigned long n_events = 100e3; 

    for (unsigned long i=0; i<n_events; i++) hist->Fill( gen_event(xmin, xmax, mrand) ); 

    new TCanvas;
    gPad->SetLogy(1); 
    hist->Draw("E"); 

    //normalize the data: 


    auto data = make_histogram_copy(hist); 

    // for the given field, scale all values so they fall in the range [-1, +1]
    auto get_maxmin = [&data](double HistoBin::*field, double& min, double& max) {
        min=+1e30; max=-1e30; 
        for (const auto& bin : data.bins) { 
            min = std::min(bin.*field, min); 
            max = std::max(bin.*field, max); 
        }
    };

    xmin=data.GetXmin(); 
    xmax=data.GetXmax(); 

    std::vector<GP::Point> all_points; 
    all_points.reserve(data.GetNbins()); 
    
    double ymin, ymax; 
    get_maxmin(&HistoBin::N, ymin, ymax); 

    ymin = std::log(ymin); 
    ymax = std::log(ymax); 

    double y_mean  = (ymax + ymin)/2.;
    double y_scale = (ymax - ymin)/2.; 

    double x_mean  = (xmax + xmin)/2.;
    double x_scale = (xmax - xmin)/2.;

    std::vector<double> pts_x, pts_y, pts_stddev; 

    for (const auto& bin : data.bins) {

        //find the bin center
        double x = (bin.xmin + bin.xmax)/2.;
        
        //scale the bin center 
        x = (x - x_mean)/x_scale; 

        // find the log of the bin contents (and the leading central moment of the uncertainty) 
        // (we're going to fit in log-space)
        double y = std::log(bin.N); 
        double stddev = 1./std::sqrt(bin.N);
        
        y = (y - y_mean)/y_scale; 
        stddev *= 1./y_scale; 
    
        all_points.emplace_back(x, y, stddev*stddev);
    
        pts_x.push_back(x);
        pts_y.push_back(y);
        pts_stddev.push_back(stddev); 
    }

    new TCanvas; 
    auto g = new TGraphErrors(pts_x.size(), pts_x.data(), pts_y.data(), nullptr, pts_stddev.data()); 
    g->Draw(); 

    std::vector<double> params{1., 1.}; 

    Kernel square_exp_kernel{
        [](double x1, double x2, const std::vector<double>& params){
            
            if (params.size() != 2) return peak_search::numbers::nan; 
            double amplitude{params[0]}, length_scale{params[1]}; 
            double arg = (x1 - x2)/length_scale; 
            return amplitude * std::exp( -0.5*arg*arg ); 
        }, params
    };

    blind_window_fitdata fitdata{
        all_points, square_exp_kernel, 8, 1
    };

    blind_window_fit(2, params.data(), nullptr, &fitdata); 
}
//_____________________________________________________________________________________________________________________
double blind_window_fit(unsigned n, const double* params, double* grad, void* _d) 
{
#ifdef DEBUG 
    static std::mutex count_mut;
    static unsigned n_calls=0;

    count_mut.lock(); ++n_calls; count_mut.unlock(); 

    Info(__func__, "<call %4u>: In body. %u params", n_calls, n); 
#endif
    const double kNaN = peak_search::numbers::nan; 

    using namespace peak_search; 

    blind_window_fitdata* fitdata = (blind_window_fitdata*)_d;
    if (!fitdata) {
        std::ostringstream oss; 
        oss << "in <" << __func__ << "> fitdata ptr is null"; 
        throw std::invalid_argument(oss.str()); 
        return kNaN; 
    }

    auto kernel = fitdata->kernel; 
    
    //give the kernel our current list of params
    kernel.params.assign(params, params+n); 

    const auto& all_points  = fitdata->data; 
    const int n_bins        = all_points.size(); 
    const int window_size   = fitdata->blind_window_size; 
    const int buffer        = fitdata->buffer; 

#ifdef DEBUG
    Info(__func__, "%i bins total, buffer=%i, window size=%i", n_bins, buffer, window_size); 
#endif

    if (n_bins < window_size + 2*buffer) {
        std::ostringstream oss; 
        oss << "in <" << __func__ << 
            "> number of bins ("<<n_bins<<") is less than window size plus buffer ("
            << window_size << " + 2*"<< buffer << ")."; 
        throw std::invalid_argument(oss.str()); 
        return kNaN; 
    }

    std::vector<GP::Point> points_sideband, points_blind, points_unblind; 

    points_sideband.reserve(n_bins - window_size); 
    points_blind  .reserve(window_size); 
    points_unblind.reserve(window_size); 

#ifdef DEBUG
    //fill out an array of points to draw
    std::vector<double> pts_x(n_bins), pts_y(n_bins), pts_err(n_bins); 
    unsigned i=0; 
    for (const auto& pt : all_points) { pts_x[i]=pt.x; pts_y[i]=pt.y; pts_err[i]=std::sqrt(pt.variance); ++i; }
#endif

    const int n_trials = n_bins - window_size - 2*buffer + 1;
#ifdef DEBUG
    Info(__func__, "starting %i trials", n_trials); 
    auto canv = new TCanvas; 
#endif

    double avg_chi2 = 0.;
    for (int i=0; i<n_trials; i++) {

        int first_bin = buffer + i; 
        int last_bin = first_bin + window_size;

#ifdef DEBUG
        Info(__func__, "trial %i/%i, bin range: [%i, %i]", i, n_trials-1, first_bin, last_bin); 
        std::vector<double> gp_x(window_size), gp_y(window_size), gp_error(window_size); 
#endif
    
        //make the sideband and blind vector of points
        points_blind  .assign( all_points.cbegin()+first_bin, all_points.cbegin()+last_bin ); 
        points_unblind.assign( all_points.cbegin()+first_bin, all_points.cbegin()+last_bin ); 
        
        points_sideband.clear(); 
        //auto dest_it = points_sideband.begin(); 
        points_sideband.insert(points_sideband.end(), all_points.cbegin(), all_points.cbegin()+first_bin); 
        points_sideband.insert(points_sideband.end(), all_points.cbegin()+last_bin, all_points.cend()); 

#ifdef DEBUG
        Info(__func__, "starting GP fit:"); 
#endif
        // Compute the Gaussian process for this point 
        GP::Compute(points_sideband, points_blind, kernel);
    
        double chi2=0.; 


        // now, evaluate the fit 
        for (int j=0; j<window_size; j++) {
            const auto& pt_blind = points_blind[j];
            const auto& pt_unblind = points_unblind[j]; 

            double arg = pt_blind.y - pt_unblind.y; 
            chi2 += (arg*arg)/(pt_blind.variance + pt_unblind.variance); 

#ifdef DEBUG
            gp_x[j]=pt_blind.x; gp_y[j]=pt_blind.y; gp_error[j]=std::sqrt(pt_blind.variance + pt_unblind.variance);
            Info(__func__, "done."); 
#endif
        }
#ifdef DEBUG
        Info(__func__, "chi2: %.4e", chi2); 

        //make the graph of all points
        auto g_pts = new TGraphErrors(n_bins, pts_x.data(), pts_y.data(), nullptr, pts_err.data()); 
        g_pts->Draw("A P Z");
        
        auto g_pred = new TGraphErrors(window_size, gp_x.data(), gp_y.data(), nullptr, gp_error.data()); 

        g_pred->SetFillColor(kGray);
        g_pred->SetLineStyle(0); 
        g_pred->Draw("SAME 3");

        auto g_pred_line = new TGraph(window_size, gp_x.data(), gp_y.data()); 
        g_pred_line->Draw("SAME L");

        canv->Modified(); 
        canv->Update(); 
        canv->SaveAs("plots/test_GP.gif+50"); 
#endif
        avg_chi2 += chi2;
    }

    return avg_chi2/((double)n_trials); 
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
//_____________________________________________________________________________________________________________________
double gen_event(double min, double max, TRandom3& gen)
{
    double x; 
    do { x = gen.Gaus()*2. - 1.; x += 0.2*std::sin(x*3.1415926536/5.); } while (x > max || x < min); 
    return x; 
}