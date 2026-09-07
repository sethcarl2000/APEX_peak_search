#ifndef peak_search_FitTest_Configuration_hpp
#define peak_search_FitTest_Configuration_hpp

#include <FitTest/ParameterList.hpp>
// stdlib
#include <cstddef> 


namespace peak_search
{
namespace FitTest
{

struct Configuration {

    /// @brief List of parameters to scan. 
    ParameterList params{}; 

    /// @brief Amount of stats for the full spectrum 
    double total_stats{75e6}; 

    /// @brief mass bin size. 
    double bin_size{0.5}; 

    /// @brief steps for a single thread to execute per task 
    unsigned long n_steps_per_task{1}; 

    /// @brief Number of threads. 0 = use all available threads
    size_t n_threads{0};
};

}
}


#endif