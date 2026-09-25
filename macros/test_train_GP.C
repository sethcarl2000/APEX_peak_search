
#include <GP.hpp>
#include <Histo1D.hpp> 
#include <make_histogram_copy.hpp>
#include <numbers.hpp>
// NLopt (nonlinear optimization lib)
#include <nlopt.hpp> 
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
// stdlib 
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

    int n_tests=200; 

    std::vector<double> param_range_amplitude{0.01, 10000.}; 
    std::vector<double> param_range_length{0.02, 2.5}; 

    auto hist_param_space = new TH2D("h_param_space", "Parameter space;log Amplitude;log Length Scale;NLL", 
        n_tests, std::log(param_range_amplitude[0]), std::log(param_range_amplitude[1]), 
        n_tests, std::log(param_range_length[0]),    std::log(param_range_length[1])
    ); 

    auto xax = hist_param_space->GetXaxis(); 
    auto yax = hist_param_space->GetYaxis(); 
 
    std::vector<double> params{std::log(1.), std::log(0.8)};

    GP::Kernel square_exp_kernel{
        [](double x1, double x2, const std::vector<double>& params){
            
            if (params.size() != 2) return peak_search::numbers::nan; 
            double amplitude{std::exp(params[0])}, length_scale{std::exp(params[1])}; 
            double arg = (x1 - x2)/length_scale; 
            return amplitude * std::exp( -0.5*arg*arg ); 
        }, params
    };

    blind_window_fitdata fitdata{
        all_points, square_exp_kernel, 8, 1
    };
    //double chi2 = blind_window_fit(2, params.data(), nullptr, &fitdata); 
    //return; 
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 

    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-5);
    minimizer->SetPrintLevel(3);    

    auto nll_computer = std::make_unique<GP::LikelihoodComputer>(all_points, square_exp_kernel); 

    double nll; 

    auto f_fcn = [&nll_computer](const double *par) { return nll_computer->ComputeNLL(par); };

    auto f_objective = ROOT::Math::Functor(f_fcn, 2); 

    minimizer->SetFunction(f_objective); 

    minimizer->SetVariable(0, "log_Amplitude", params[0], 1e-4); 
    minimizer->SetVariable(1, "log_Length",    params[1], 1e-5); 
    
    auto fit_result = minimizer->Minimize(); 

    if (!fit_result) {
        Error(__func__, "Something went wrong with the fit result!");   
        return; 
    }

    //copy parameter data 
    if (auto pars = minimizer->X(); pars!=nullptr) 
        params.assign(pars, pars+2); 

    nll = minimizer->MinValue(); 

    std::printf("optimal params: amplitude=%.1f, length=%.4f      NLL=%.3e\n", 
        std::exp(params[0]), std::exp(params[1]), nll
    ); 

    std::vector<double> cov_mat_elems(params.size()*params.size(), 0.);
    
    //get the coviariance matrix 
    Eigen::MatrixXd cov(params.size(), params.size()); 
    minimizer->GetCovMatrix(cov.data()); 

    //get eigenvectors / eigenvalues 
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(cov); 

    if (eigen_solver.info() != Eigen::Success) {
        Error(__func__, "Something went wrong finding eigenvectors of covariance matrix!"); 
        return; 
    }
    Eigen::MatrixXd e_vecs = eigen_solver.eigenvectors(); 
    Eigen::VectorXd e_vals = eigen_solver.eigenvalues(); 

    std::cout << "cov. matrix: \n"; 
    for (int i=0; i<params.size(); i++) {
        for (int j=0; j<params.size(); j++) {
            std::cout << cov(i,j) << " "; 
        }
        std::cout << "\n"; 
    }

    std::cout << "Eigenvalues / vectors: \n"; 
    for (int i=0; i<params.size(); i++) {
        std::cout << " " << e_vals(i) << "      { "; 
        for (int j=0; j<params.size(); j++) std::cout << e_vecs(i,j) <<  " "; 
        std::cout << "}\n"; 
    }

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
        g_pts->SetMarkerStyle(kPlus); 
        g_pts->Draw("A P Z");
        
        auto g_pred = new TGraphErrors(window_size, gp_x.data(), gp_y.data(), nullptr, gp_error.data()); 

        g_pred->SetFillColorAlpha(kGray, 0.4); 
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
double marginal_nll(unsigned n, const double* params, double* grad, void* data)
{
    using namespace peak_search; 

    auto computer = (GP::LikelihoodComputer*)data; 

    if (!computer) {
        throw std::invalid_argument("in <marginal_nll>: computer ptr is null"); 
        return peak_search::numbers::nan; 
    }
    return computer->ComputeNLL(params); 
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