
#include <GP.hpp>
// stdlib
#include <cmath> 

namespace peak_search
{
namespace GP
{

//__________________________________________________________________________________________________________________
LikelihoodComputer::LikelihoodComputer(const std::vector<Point>& points, Kernel kernel)
        : fPoints{points}, fKernel{kernel}, n_params{(unsigned)kernel.params.size()}, n_pts{(unsigned)points.size()} 
{
    if (points.empty()) {
        throw std::logic_error("no points given to fit!"); 
        return; 
    }
    
    if (n_params < 0) {
        throw std::logic_error("in <LikelihoodComputer::ctor> this kernel does not have any params!"); 
        return; 
    }

    //allocate the covariance matrix
    fCovMatrix = Eigen::MatrixXd(n_pts, n_pts); 
    fLLT = Eigen::LLT<Eigen::MatrixXd>(n_pts); 

    // vector of y-values
    fY = Eigen::VectorXd(n_pts); 
    for (unsigned i=0; i<n_pts; i++) fY(i) = fPoints[i].y; 

    fPrefactor = 0.5 * ((double)n_pts) * std::log(2.*3.1415926536); 
}
//__________________________________________________________________________________________________________________    
double LikelihoodComputer::ComputeNLL(const double* params)
{

    // copy this set of parameters over to our kernel function
    fKernel.params.assign(params, params+n_params); 

    //now, compute the covariance matrix
    for (unsigned i=0; i<n_pts; i++) { 
     
        double xi = fPoints[i].x; 

        // add the diagonal variance (and the error of this datapoint)
        fCovMatrix(i,i) = fKernel(xi, xi) + fPoints[i].variance; 

        for (unsigned j=i+1; j<n_pts; j++) {
     
            double xj = fPoints[j].x; 
            fCovMatrix(i,j) = fKernel(xi, xj); 
        }
    }

    // symmetrize the matrix
    for (unsigned i=1; i<n_pts; i++) 
        for (unsigned j=0; j<i; j++) fCovMatrix(i,j) = fCovMatrix(j,i); 

    //now, solve the system
    fLLT.compute(fCovMatrix); 

    double negative_log_likelihood =0.; 

    negative_log_likelihood += 0.5 * fY.dot( fLLT.solve(fY) ); 

    // 1/2 log det|CovMatrix| = 1/2 log( det(|L|)^2 ) = sum_i log(L_ii)
    double log_det = fLLT.matrixLLT().diagonal().array().log().sum(); 

    negative_log_likelihood += log_det; 

    return negative_log_likelihood + fPrefactor; 
}

}
}