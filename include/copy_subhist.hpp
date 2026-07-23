#ifndef peak_search_copy_subhist_hpp
#define peak_search_copy_subhist_hpp

#include <bininfo.hpp>
#include <TH1D.h>

namespace peak_search
{

/// @brief Copies sub-range of histogram. if 'xmin' or 'xmax' slices a particular bin, that whole bin is included. 
/// @param data histogram data to copy
/// @param xmin start of x-range to copy
/// @param xmax end of x-range to copy
/// @return Copy of histogram over given sub-range
histo_1D_t copy_subhist(TH1D* hist, double xmin, double xmax);

};

#endif