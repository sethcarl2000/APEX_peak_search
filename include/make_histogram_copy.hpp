#include <Histo1D.hpp>
#include <numbers.hpp>
// ROOT headers
#include <TH1D.h> 

namespace peak_search {

/// @brief Make a copy of a TH1D (as a Histo1D)
/// @param hist_source source histogram to  get the data from
/// @param xmin min x-value to include (if it cuts a bin in half, the whole bin is included). if NaN is passed, or outside source hist range, lower bound is source hist lower bound. 
/// @param xmax max x-value to include (if it cuts a bin in half, the whole bin is included). if NaN is passed, or outside source hist range, upper bound is source hist upper bound. 
/// @return 
Histo1D make_histogram_copy(TH1D* hist_source, double xmin=numbers::nan, double xmax=numbers::nan);

};