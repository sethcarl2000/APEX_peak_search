#ifndef peak_search_GP_Compute_hpp
#define peak_search_GP_Compute_hpp

#include <GP/Point.hpp>
#include <Fcn1D/Fcn1D.hpp>
// stdlib headers
#include <vector> 

namespace peak_search
{
namespace GP
{

/// @brief Given a set of input points, output points, and kernel, compute the GP mean and variance estimate for each output point. 
/// @param inputs input points, with given uncertainties
/// @param outputs output points; for each x-value, we will estimate mean & variance  
/// @param kernel kernel to use. input is the distance between the x-value of two points, output is the covariance between their outputs: k(x_i, x_j) = Cov(y_i, y_j)
void Compute(const std::vector<Point>& inputs, std::vector<Point>& outputs, const Fcn1D& kernel); 

}
}


#endif