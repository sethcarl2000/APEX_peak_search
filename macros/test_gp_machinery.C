
#include "cov_matrix_to_ellipse.h"
#include <GP.hpp>
#include <Histo1D.hpp> 
#include <make_histogram_copy.hpp>
#include <numbers.hpp>
// analysis_utils lib
#include <analysis_utils/task.hpp>
#include <analysis_utils/ROOT.hpp>
#include <analysis_utils/thread_pool.hpp>
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
#include <Math/ProbFuncMathCore.h>
#include <TLine.h> 
#include <TError.h> 
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TStopwatch.h> 
// stdlib 
#include <thread> 
#include <vector> 
#include <cmath> 
#include <functional> 
#include <mutex> 
#include <memory> 


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

    const std::size_t n_pts = hist->GetXaxis()->GetNbins(); 

    //normalize the data: 
    std::vector<GP::Point> log_data; log_data.reserve(n_pts); 

    auto hist_data = make_histogram_copy(hist); 

    for (const auto& bin : make_histogram_copy(hist).bins) {
        log_data.emplace_back(
            0.5*(bin.xmin + bin.xmax), //x-position
            std::log10(bin.N), //value
            1./( bin.N * std::log(10.)*std::log(10) ) //variance
        ); 
    }   

    GP::TrainParameters params; 

    params.draw_fit = false; 

    // log-data 
    params.data = log_data; 

    std::vector<GP::HyperParameter> gp_params{
        {1., {-2., +6.}, "log_rms_slope", 0.5},
        {0., {-6., +6.}, "log_A",         0.5}
    };  


    /// Param definitions: 
    /// par[0] = log RMS slope log(A / sigma^2)
    /// par[1] = log characteristic length scale log(sigma)
    std::vector<double> pars{1., 0.}; 
    params.kernel = GP::Kernel{ 
        [](double x1, double x2, const std::vector<double>& params){
            
            if (params.size() != 2) return peak_search::numbers::nan; 

            double sigma = std::exp( 0.5*params[1] - params[0] );
            double A     = std::exp( params[1] ); 

            double arg = (x1 - x2)/sigma; 
            return A * std::exp( -0.5*arg*arg ); 
        }, {1., 0.}
    };
    
    params.params = gp_params; 
    params.minuit_printlevel = 2; 

    params.draw_fit = true; 

    auto computer = std::make_unique<GP::Computer>(); 

    GP::Train(params, computer.get()); 

    const bool fit_success = static_cast<bool>(params); 

    if (!fit_success) {
        Error(__func__, "Something went wrong with the GP training. what(): %s", params.status.c_str()); 
    }
        
    std::vector<double> pts_x(n_pts), pts_y(n_pts), pts_err(n_pts), pts_predict(n_pts), pts_residual(n_pts);

    for (std::size_t i=0; i<n_pts; i++) {

        pts_x[i] = log_data[i].x; 
        pts_y[i] = log_data[i].y; 
        pts_err[i] = std::sqrt( log_data[i].variance ); 
    }

    auto log_data_predict = log_data; 

    computer->Compute(log_data, log_data_predict); 

    for (std::size_t i=0; i<n_pts; i++) {
        pts_x[i] = log_data[i].x;
        pts_y[i] = log_data[i].y; 
        pts_predict[i] = log_data_predict[i].y; 
        
        pts_residual[i] = (log_data[i].y - log_data_predict[i].y)/std::sqrt(log_data[i].variance); 
    }

    TCanvas *canv;
    TGraph *g;
    TGraphErrors *ge; 

    canv = new TCanvas; 
    g = new TGraph(n_pts, pts_x.data(), pts_y.data()); 
    
    g->SetTitle("Histogram and GP compatibility;mass;log_{10} bin count"); 
    
    g->SetMarkerStyle(kPlus); 
    g->SetMarkerSize(0.5); 
    g->Draw("AP"); 
    
    if (fit_success) {
        g = new TGraph(n_pts, pts_x.data(), pts_predict.data()); 
        g->SetLineStyle(kSolid); 
        g->SetLineColor(kRed);
        g->Draw("SAME"); 

        auto std_err = std::vector<double>(n_pts, 1.); 

        canv = new TCanvas; 
       
        g = new TGraph(n_pts, pts_x.data(), pts_residual.data()); 
        g->SetTitle("Per-bin residuals from GP fit;mass;(x_{i} - <x_{i}>) / #sigma_{i}"); 
        g->SetMarkerStyle(kFullCircle); 
        g->SetFillColor(kBlack);
        g->SetMarkerSize(0.5); 
        g->Draw("APL"); 

        auto l = new TLine; 
        l->SetLineColor(kRed); 
        l->DrawLine(pts_x.front(),0., pts_x.back(), 0.); 
        
        l->SetLineStyle(kDashed); 
        l->DrawLine(pts_x.front(),+1., pts_x.back(),+1); 
        l->DrawLine(pts_x.front(),-1., pts_x.back(),-1);

        auto hist_residual = new TH1D("h_residual", "P-value distribution of bin-residuals (GP Fit);p(residual);", 25, 0., 1.); 
        for (auto z : pts_residual) {
            hist_residual->Fill( ROOT::Math::normal_cdf(z, 1., 0.) ); 
        }
        canv = new TCanvas; 
        hist_residual->SetMinimum(0.); 
        hist_residual->SetMaximum( hist_residual->GetMaximum()*3. ); 
        hist_residual->Draw("HIST E"); 
    }    
    return; 
    
    Eigen::MatrixXd cov_matrix = params.covariance_matrix; 

    double L1, L2, theta; 
    if (fit_success) {
        cov_matrix_to_ellipse(
            cov_matrix(0,0), cov_matrix(1,1), cov_matrix(0,1), 
            L1, L2, theta
        );
    }

    double fit_x{params.kernel.params[0]}, fit_y{params.kernel.params[1]}; 

    // scan 2d parameter space 

    auto all_points = GP::Normalize_points(params.data); 
    auto nll_computer = std::make_unique<GP::LikelihoodComputer>(all_points, params.kernel); 

    int n_bins = 200; 

    auto hist_param_space = new TH2D(
        "h_data", 
        "Parameter space for k(|x|) = A e^{-|x|^{2}/2#sigma^{2} };log[ #sqrt{ <f'(x)^{2}> } ] = log[ #sqrt{A} / #sigma ];log[ A ];-log L",  
        n_bins, -2, 6., 
        n_bins, std::log(0.001), std::log(800.)
    ); 
    double dx = hist_param_space->GetXaxis()->GetBinWidth(1); 
    double dy = hist_param_space->GetYaxis()->GetBinWidth(1); 

    double max_nll = 2.5e3; 

    auto xax = hist_param_space->GetXaxis();
    auto yax = hist_param_space->GetYaxis(); 

    struct nll_pt { double x, y, val; };

    using namespace analysis_utils;    
    auto plist = task::param_list<nll_pt>(); 

    plist.add(&nll_pt::x, xax->GetBinCenter(1), xax->GetBinCenter(n_bins), n_bins); 
    plist.add(&nll_pt::y, yax->GetBinCenter(1), yax->GetBinCenter(n_bins), n_bins); 

    task::config_t cfg; 

    cfg.chunk_size = 32; 
    cfg.verbosity = 1;
    //cfg.n_threads = 32; 

    auto t_hist = root::threadlocal_hist(hist_param_space); 

    thread_pool pool; 

    auto t_cptr = pool.MakeThreadLocalObj(*nll_computer.get()); 
    auto t_pars = pool.MakeThreadLocalObj(pars); 

    TStopwatch timer; 
    task::scan_params(plist, [&t_cptr,&t_pars,max_nll,&t_hist](const nll_pt *__restrict point, size_t t){

        //
        auto& mpars = t_pars(t);
        auto& my_cptr = t_cptr(t); 
        if (mpars.empty()) return; 
        mpars[0] = point->x; 
        mpars[1] = point->y; 
        double nll = my_cptr.ComputeNLL(mpars.data());   
        if (!(peak_search::numbers::is_nan(nll) || nll > max_nll)) { 
            t_hist(t)->Fill( point->x, point->y, nll ); 
        }
        
    }, cfg); 

    auto real_time = timer.RealTime(); 
    std::printf("time taken for %zi steps: %.4f (%.4f ms / step)\n", 
        plist.get_n_steps(), 
        real_time, 
        (real_time*1e3)/((double)plist.get_n_steps())
    ); 

    //combine all histograms into one 
    hist_param_space = t_hist.aggregate(); 

    canv = new TCanvas; 
    canv->SetRightMargin(0.15); 
    hist_param_space->Draw("colz"); 

    if (fit_success) {
        std::printf(
            "center: (%+10.4e, %+10.4e)\n"
            "eigenvalues / eigenvectors:\n" 
            "   %+10.4e     [%+10.4e, %+10.4e]\n"
            "   %+10.4e     [%+10.4e, %+10.4e]\n"
            "theta: %+.1f deg\n",
            fit_x, fit_y, 
            L1, std::cos(theta),  std::sin(theta), 
            L2, std::sin(-theta), std::cos(theta),
            theta*(180./3.14159265)
        );
        auto ellipse = new TEllipse(fit_x, fit_y, L1, L2, 0., 360., (180./3.14159265)*theta); 
        ellipse->SetFillStyle(0);
        ellipse->SetLineColor(kRed); 
        ellipse->Draw("SAME"); 

        ellipse = new TEllipse(fit_x, fit_y, 2.*L1, 2.*L2, 0., 360., (180./3.14159265)*theta); 
        ellipse->SetFillStyle(0);
        ellipse->SetLineColor(kRed);
        ellipse->SetLineStyle(kDashed);  
        ellipse->Draw("SAME"); 

        TLine* line = new TLine; 
        line->SetLineColor(kRed); 
        line->SetLineStyle(kDashed); 

        line->DrawLine(fit_x,yax->GetXmin(), fit_x,yax->GetXmax()); 
        line->DrawLine(xax->GetXmin(),fit_y, xax->GetXmax(),fit_y); 
    }
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
template<typename F, int Nthreads> void parallel_for(std::size_t n_tasks, const F& fcn, std::size_t chunk_size, std::size_t n_threads)
{
    //keeps track of tasks
    alignas(64) std::atomic<std::size_t> task_id{0}; 

    // if there are no tasks, then quit
    if (n_tasks<1) return; 

    //get the number of threads
    if (n_threads<1) n_threads = std::thread::hardware_concurrency(); 

    std::array<std::thread, Nthreads> threads; threads.reserve(n_threads); 

    if (chunk_size<1) chunk_size = std::max<std::size_t>( 1, n_tasks / (100 * n_threads)); 

    for (size_t t=0; t<n_threads; t++) {

        threads[t] = [&fcn, &task_id, n_tasks, chunk_size, n_threads, t]{

            while (1) {
                std::size_t start = task_id.fetch_add(chunk_size, std::memory_order_relaxed); 
                std::size_t end   = std::min<std::size_t>( n_tasks, start + chunk_size ); 

                //quit if all tasks have already been scheduled
                if (start >= n_tasks) break;    

                for (std::size_t index=start; index<end; index++) fcn(index, t); 
            }; 
        }; 
    }

    for (auto& thread : threads) thread.join(); 
}
//_____________________________________________________________________________________________________________________
double gen_event(double min, double max, TRandom3& gen)
{
    double x; 
    do { x = gen.Gaus()*2. - 1.; x += 0.2*std::sin(x*3.1415926536/5.); } while (x > max || x < min); 
    return x; 
}