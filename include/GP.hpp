#ifndef peak_search_GP_hpp
#define peak_search_GP_hpp

#include <Fcn1D/Fcn1D.hpp>
#include <numbers.hpp>
#include <Histo1D.hpp>
// Eigen 
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
// stdlib
#include <vector> 
#include <functional> 
#include <memory> 
#include <string> 

namespace peak_search
{
namespace GP
{

/// @brief A single point for regression in a gaussian process
struct Point { double x, y, variance; };

/// @brief Kernel function 
struct Kernel {

    std::function<double(double,double,const std::vector<double>&)> fcn{
        // this is just a dummy placeholder default kernel function, to make this struct default-constructible 
        [](double,double,const std::vector<double>&){ return peak_search::numbers::nan; }
    }; 

    std::vector<double> params{}; 

    double operator()(double x1, double x2) const { return fcn(x1,x2,params); }
}; 
//check that this struct is default constructible 
static_assert(std::is_default_constructible_v<Kernel>);


/// @brief Given a set of input points, output points, and kernel, compute the GP mean and variance estimate for each output point. 
/// @param inputs input points, with given uncertainties
/// @param outputs output points; for each x-value, we will estimate mean & variance  
/// @param kernel kernel to use. input is the distance between the x-value of two points, output is the covariance between their outputs: k(x_i, x_j) = Cov(y_i, y_j)
void Compute(const std::vector<Point>& inputs, std::vector<Point>& outputs, const Fcn1D& kernel); 


/// @brief Given a set of input points, output points, and kernel, compute the GP mean and variance estimate for each output point. 
/// @param inputs input points, with given uncertainties
/// @param outputs output points; for each x-value, we will estimate mean & variance  
/// @param kernel kernel to use. input is the distance between the x-value of two points, output is the covariance between their outputs: k(x_i, x_j) = Cov(y_i, y_j)
void Compute(const std::vector<Point>& inputs, std::vector<Point>& outputs, const Kernel& kernel); 

/// @brief Given a set of observations (with given gaussian errors), comptue the marginal likelihood 
/// @param points set of observations to compute the likelihood w/r/t   
/// @param kernel kernel to use. input is the distance between the x-value of two points, output is the covariance between their outputs: k(x_i, x_j) = Cov(y_i, y_j)
double MarginalLikelihood(const std::vector<Point>& points, const Kernel& kernel);

/// @brief Computes marginal likelihood w/r/t a certain set of points 
class LikelihoodComputer {
private: 
    std::vector<Point> fPoints{}; 
    Kernel fKernel{}; 
    unsigned n_params, n_pts; 

    // matrices needed for covariance / lin system solving. 
    Eigen::MatrixXd fCovMatrix; 
    Eigen::LLT<Eigen::MatrixXd> fLLT; 
    Eigen::VectorXd fY; 

    // prefactor for marginal NLL 
    double fPrefactor; 

public: 
    LikelihoodComputer(const std::vector<Point>& points, Kernel kernel); 

    /// @brief Compute marginal likelihood w/r/t a given set of points 
    /// @param pars ptr to parameter array 
    /// @return (Negative) Marginal log-likelihood of points w/r/t current kernel parameters 
    double ComputeNLL(const double* pars); 
};
static_assert(std::is_copy_constructible_v<LikelihoodComputer>); 

/// @brief Given a histogram, return a normalized set of GP::Point-s (with y -> ln(y)), so that all x-values and y-values span the range [-1,1]
/// @param data input poitns to normlaize
/// @return output datapoints, normalized w/r/t x and y s.t. x \in [-1,1], ln(y) \in [-1,1]
std::vector<Point> Normalize_points(const Histo1D& data); 

struct TrainResult {

    Eigen::MatrixXd param_covariance; 
    std::vector<double> params; 
    bool is_success{false}; 

    operator bool() const { return is_success; }
};

template<typename T> struct Bound {
    T min, max; 
    T clamp(const T& val) const { return (val<min ? min : (val>max ? max : val)); }
}; 

struct TrainParameters {

    std::string status{"null"}; 
    Kernel kernel; 
    Eigen::MatrixXd covariance_matrix; 
    Histo1D data; 
    std::vector<Bound<double>> params_bounds{}; 
    bool draw_fit{false}; 

    operator bool() const { return status == "Success"; }
};

//class that actually does GP on a set of points
class Computer {
private: 

    friend void Train(TrainParameters&, Computer*); 

    Kernel fKernel{}; 

    double x_scale, x_mean; 
    double y_scale, y_mean; 

public: 

    Computer() = default; 

    /// @brief Compte GP with given kernel 
    /// @param inputs input (sideband) points
    /// @param outputs output (blind) points
    void Compute(std::vector<Point> inputs, std::vector<Point>& outputs) const; 

};

/// @brief Train the hyperparameters of the given GP 
/// @param kernel kernel to use 
void Train(TrainParameters& params, Computer* computer=nullptr); 

}
}




#endif