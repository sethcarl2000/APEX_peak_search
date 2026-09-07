#ifndef peak_search_FitTest_Outputs_hpp
#define peak_search_FitTest_Outputs_hpp

#include <FitTest/ParameterList.hpp>
// ROOT
#include <TObject.h>
// stdlib
#include <cstddef> 
#include <vector> 


namespace peak_search
{
namespace FitTest
{

class Outputs {
private: 

    std::vector<TObject*> fOutputs{}; 

public: 

    Outputs() = default;

    size_t Add(TObject* obj) { fOutputs.emplace_back(obj); return fOutputs.size()-1; } 

    std::vector<TObject*>& GetPtrs() { return fOutputs; }
};

}
}


#endif