
#include <GP.hpp>
#include <Histo1D.hpp>
// ROOT
#include <TCanvas.h>
#include <TGraph.h> 
#include <TGraphErrors.h>
#include <TH2D.h> 
#include <TString.h> 
#include <TError.h> 
// ROOT/Minuit 
#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
// stdlib
#include <stdexcept> 


namespace peak_search
{
namespace GP
{

/// @brief Train the hyperparameters of the given GP 
/// @param data data to train on
/// @param kernel kernel to use 
void Train(TrainParameters& params, Computer* computer)
{
    auto& data = params.data; 

    // for the given field, scale all values so they fall in the range [-1, +1]
    auto get_maxmin = [&data](const std::vector<Point>& pts, double Point::*field, double& min, double& max) {
        min=+1e30; max=-1e30; 
        for (const auto& pt : pts) { 
            min = std::min(pt.*field, min); 
            max = std::max(pt.*field, max); 
        }
    };

    bool make_computer = (computer != nullptr); 

    // the actual kernel, whose hyperparamters we will seek to optimize
    auto& kernel = params.kernel;   
    if (make_computer) computer->fKernel = kernel; 

    double xmin=data.GetXmin(); 
    double xmax=data.GetXmax(); 

    auto all_points = Normalize_points(params.data); 
    
    double ymin, ymax; 
    get_maxmin(all_points, &Point::y, ymin, ymax); 

    ymin = std::log(ymin); 
    ymax = std::log(ymax); 

    double y_mean  = (ymax + ymin)/2.;
    double y_scale = (ymax - ymin)/2.; 

    double x_mean  = (xmax + xmin)/2.;
    double x_scale = (xmax - xmin)/2.;

    // if true, we will draw the fit. 
    const bool do_draw = params.draw_fit; 

    const size_t n_params = kernel.params.size(); 

    if (make_computer) {

        computer->x_mean  = x_mean; 
        computer->x_scale = x_scale; 
        computer->y_mean  = y_mean; 
        computer->y_scale = y_scale; 
    }

    const size_t n_points = all_points.size(); 
    std::vector<double> pts_x, pts_y, pts_stddev; 
    TCanvas* canv; 
    if (do_draw) {
        pts_x       .reserve(n_points); 
        pts_y       .reserve(n_points); 
        pts_stddev  .reserve(n_points);
        for (const auto& pt : all_points) {
            
            if (do_draw) {
                pts_x.push_back(pt.x);
                pts_y.push_back(pt.y);
                pts_stddev.push_back(std::sqrt(pt.variance)); 
            }
        }
        canv = new TCanvas; 
        auto g = new TGraphErrors(pts_x.size(), pts_x.data(), pts_y.data(), nullptr, pts_stddev.data()); 
        g->Draw(); 
    }
    
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 

    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-5);
    minimizer->SetPrintLevel(1);    

    auto nll_computer = std::make_unique<GP::LikelihoodComputer>(all_points, kernel); 

    double nll; 

    auto f_fcn = [&nll_computer](const double *par) { return nll_computer->ComputeNLL(par); };

    auto f_objective = ROOT::Math::Functor(f_fcn, 2); 

    minimizer->SetFunction(f_objective); 

    auto& pars = kernel.params; 

    unsigned i_par=0; 
    for (const auto par : pars) {
        minimizer->SetVariable(i_par, Form("par_%u",i_par), par, 1e-4); 
        ++i_par; 
    }

    Info(__func__, "Starting minimization..."); 
    auto fit_result = minimizer->Minimize();

    if (!fit_result) {

        Error(__func__, "Something went wrong with the fit result!");   
        params.status = "Minimization failed"; 
        return; 
    }
    Info(__func__, "Done."); 

    //copy parameter data 
    if (auto ptr = minimizer->X(); ptr!=nullptr) {

        kernel.params.assign(ptr, ptr+n_params);        
 
    } else {

        Error(__func__, "Array of best-fit parameters is null");   
        params.status = "Array of best-fit parameters is null"; 
        return; 
    }

    nll = minimizer->MinValue(); 

    //get the coviariance matrix 
    params.covariance_matrix = Eigen::MatrixXd(n_params, n_params); 
    minimizer->GetCovMatrix(params.covariance_matrix.data()); 

    if (make_computer) {
        computer->fKernel = kernel; 
    }

    params.status = "Success"; 
}


}
}
