#ifndef peak_search_FitTest_Run_hpp
#define peak_search_FitTest_Run_hpp

#include <FitTest/Configuration.hpp>
#include <FitTest/Outputs.hpp>
#include <FitTest/Function.hpp>

namespace peak_search
{
namespace FitTest
{

/// @brief Execute the specified number of trials, given by user-provided configuration. 
/// @param n_trials number of trials to execute  
/// @param config configuration of one trial
/// @param outputs list of outputs to make available to the user function
/// @param run_verbosity Verbosity level for run executon. 0=quiet 1=print active status bar, 2=print statement at each new task submission  
/// @param fcn user-provided function to run 
void Run(size_t, Configuration, Outputs, Function, int run_verbosity=1); 

}
}


#endif