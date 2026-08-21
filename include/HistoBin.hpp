#ifndef peak_search_HistoBin_hpp
#define peak_search_HistoBin_hpp

#include <type_traits> 

namespace peak_search
{
    
struct HistoBin { double xmin, xmax, N; };

//make sure the compiler is going to generate copy / move constructors for free 
static_assert(std::is_trivially_copy_constructible_v<HistoBin>);
static_assert(std::is_trivially_move_constructible_v<HistoBin>);

};

#endif