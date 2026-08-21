#ifndef peak_search_Histo1D_hpp
#define peak_search_Histo1D_hpp



#include <HistoBin.hpp>
// stdlib headers
#include <vector> 
#include <limits> 

namespace peak_search 
{

//very simple histogram-like container 
struct Histo1D {

    std::vector<HistoBin> bins; 
    
    double GetXmax() const { return bins.empty() ? std::numeric_limits<double>::quiet_NaN() : bins.back().xmax; }
    double GetXmin() const { return bins.empty() ? std::numeric_limits<double>::quiet_NaN() : bins.front().xmin; }

    size_t GetNbins() const { return bins.size(); }
};

//check to make sure it's default constructable
static_assert(std::is_default_constructible_v<Histo1D>); 
//check to make sure we are move constructable / assignable (with no expcetion guranteed)
static_assert(std::is_nothrow_move_constructible_v<Histo1D>); 
static_assert(std::is_nothrow_move_assignable_v<Histo1D>);
//check to make sure we are copy constructable / assignable
static_assert(std::is_copy_constructible_v<Histo1D>); 
static_assert(std::is_copy_assignable_v<Histo1D>);

}; 


#endif