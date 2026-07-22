#ifndef peak_search_fit_parameter_t_hpp
#define peak_search_fit_parameter_t_hpp

#include <vector> 
#include <string> 

struct fit_parameter_t {
    double val; 
    std::string name{"param"}; 
    bool is_fixed{false}; 
};

#endif