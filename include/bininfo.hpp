#ifndef bininfo_hpp
#define bininfo_hpp

#include <TH1D.h>
#include <TH2D.h> 

#include <vector>

namespace peak_search 
{

struct bininfo_1D_t { double x, N; };

struct bininfo_2D_t { double x,y, N; };

struct histo_1D_t {
    std::vector<bininfo_1D_t> bins; 
    double xmin, xmax; 
};

struct histo_2D_t {
    std::vector<bininfo_2D_t> bins; 
    double xmin, xmax;
    double ymin, ymax;  
};

/// @brief Copy contents of histogram bins into vector 
histo_1D_t copy_1D_hist(TH1D* hist);

/// @brief Copy contents of histogram bins into vector 
histo_2D_t copy_2D_hist(TH2D* hist);

};

#endif