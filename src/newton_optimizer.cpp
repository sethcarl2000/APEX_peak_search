
#include <newton_optimizer.hpp>
#include <numbers.hpp>
#include <gauss_integrate.hpp> 

//Eigen
#include <eigen3/Eigen/Core> 
#include <eigen3/Eigen/Dense>
//ROOT
#include <Math/ProbFuncMathCore.h>
//stdlib
#include <sstream> 
#include <stdexcept> 
#include <string> 
#include <cmath> 

namespace peak_search 
{

double newton_optimizer(const histo_1D_t& data, Fcn1D& fcn, std::vector<fit_parameter_t>& params, int max_iterations)
{
    //check to make sure that the fcn and params match
    if (fcn.GetDoF() != (int)params.size()) {
        std::ostringstream oss; 
        oss << "in <"<<__func__<<">: number of fcn parameters ("<<fcn.GetDoF()<<") "
            "does not match size of parameter list ("<<params.size()<<")"; 
        
        throw std::logic_error(oss.str()); 
        return numbers::nan; 
    }

    using Eigen::MatrixXd, Eigen::VectorXd;
    
    const double dx = (data.xmax - data.xmin)/((double)data.bins.size()); 

    //set the parameters
    fcn.SetParams(params); 

    size_t n_mutable=0; 
    std::vector<size_t> ind; ind.reserve(params.size()); 
    for (size_t i=0; i<params.size(); i++) { if (!params[i].is_fixed) { ind.push_back(i); } }

    if (n_mutable < 1) return 0.; 

    double eta=0.; 

    for (int it=0; it<max_iterations; it++) {

        //loop over all bins. compute 'eta' 
        VectorXd dEta = VectorXd::Zero(ind.size()); 
        MatrixXd J    = MatrixXd::Zero(ind.size(), ind.size()); 
        
        double chi2 =0.; 
        eta = 0.; 

        std::vector<double> dL_dTi(ind.size(), 0.);             
        std::vector<double> dL_dTi_dTj(ind.size()*ind.size(), 0.); 

        for (const auto& bin : data.bins) {
            
            double x0 = bin.x - dx/2.;
            double x1 = bin.x + dx/2.; 

            double n_i = bin.N; 
            double lambda_i = gauss_integrate(fcn, x0,x1); 

            double arg = (lambda_i - n_i)/lambda_i; 

            chi2 += arg*arg*lambda_i; 
            
            eta += lambda_i - n_i*std::log(lambda_i);

            //get the derivatives of the expectation values for each bin 
            for (int i=0; i<ind.size(); i++) {

                dL_dTi[i] 
                    = gauss_integrate([&ind,i,&fcn](double x){ return fcn.Di(x, ind[i]); }, x0,x1); 

                for (int j=i; j<ind.size(); j++) {

                    dL_dTi_dTj[i*n_mutable + j] 
                        = gauss_integrate([&ind,i,j,&fcn](double x){ return fcn.Di_Dj(x, ind[i],ind[j]); }, x0,x1); 
                }
            }
            
            for (int i=0; i<ind.size(); i++) {

                dEta(i) += arg*dL_dTi[i]; 

                for (int j=i; j<n_mutable; j++) J(i,j) 
                    += dL_dTi[i]*dL_dTi[j]*(1. - arg)/lambda_i
                    +  arg*dL_dTi_dTj[i*n_mutable + j]; 
            }

        }// for (const auto& bin : data.bins)

        for (int i=1; i<n_mutable; i++) { for (int j=0; j<i; j++) J(j,i) = J(i,j); }

        VectorXd dX = J.llt().solve(dEta); 

        //add this result to the overall value
        for (size_t i=0; i<ind.size(); i++) {
            params[ind[i]].val += -dX(i);  
        }
        std::printf("<%s>: it %2i/%i, eta = %.4e, chi^2 = %.4e p(chi^2) = %.4e\n", __func__, 
            it,max_iterations, 
            eta, 
            chi2, 
            ROOT::Math::chisquared_cdf(chi2, data.bins.size())
        );

        fcn.SetParams(params); 
    }
    return eta; 
} 

};