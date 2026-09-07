#ifndef peak_search_FitTest_Function_hpp
#define peak_search_FitTest_Function_hpp

#include <SignalFit.hpp>
// stdlib 
#include <functional> 

namespace peak_search
{
namespace FitTest
{

class ThreadManager; 

using Function = std::function<void(ThreadManager*)>; 

}
}

#endif