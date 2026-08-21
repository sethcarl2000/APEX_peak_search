#ifndef peak_search_FitTestFunction_hpp
#define peak_search_FitTestFunction_hpp

#include <SignalFit.hpp>
// stdlib 
#include <functional> 

namespace peak_search
{

class FitTestThreadManager; 

using FitTestFunction = std::function<SignalFit(FitTestThreadManager*)>; 

};

#endif