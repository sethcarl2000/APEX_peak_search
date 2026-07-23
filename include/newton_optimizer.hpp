#ifndef peak_search_newton_optimizer_hpp
#define peak_search_newton_optimizer_hpp

#include <Fcn1D/Fcn1D.hpp>
#include <bininfo.hpp>
#include <fit_parameter.hpp> 

#include <vector> 

namespace peak_search
{

/// @brief Fit a function to a 1D histogram using a newtons-method optimizer. Parameters must be already 'close' for it to work
/// @param data data to fit
/// @param fcn Fcn1D fcn to use
/// @param params parameters of the function to use
/// @param max_iterations maximum number of iterations to execute
/// @return returns the 'eta' (NLL without combinatoric factor)
double newton_optimizer(const histo_1D_t& data, Fcn1D& fcn, std::vector<fit_parameter_t>& params, int max_iterations=8); 

};

#endif