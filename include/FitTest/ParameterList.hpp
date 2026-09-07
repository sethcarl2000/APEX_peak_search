#ifndef peak_search_FitTest_ParameterList_hpp
#define peak_search_FitTest_ParameterList_hpp

#include <FitTest/Parameter.hpp>
// stdlib
#include <vector>
#include <cstddef> 

namespace peak_search
{
namespace FitTest
{


class ParameterList {
private: 

    std::vector<Parameter> fParams{}; 
    std::vector<unsigned long> fStepSizes{}; 

    void ComputeStepSizes();

public: 

    ParameterList(const std::vector<Parameter>& params={}) : fParams{params}, fStepSizes{} {
        ComputeStepSizes(); 
    }

    size_t Append(const Parameter& par) {
        fParams.emplace_back(par); 
        ComputeStepSizes(); 
        return fParams.size()-1; 
    }
    
    size_t Append(unsigned long n_steps, double min, double max) {
        fParams.emplace_back(n_steps, min, max); 
        ComputeStepSizes(); 
        return fParams.size()-1; 
    }

    unsigned long GetNSteps() const { return (fParams.empty()) ? 0 : fStepSizes.back() * fParams.back().get_n_steps(); } 

    std::vector<double> GetParamList(unsigned long step) const; 

};


}
}


#endif