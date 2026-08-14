#ifndef peak_search_generate_toy_events_hpp
#define peak_search_generate_toy_events_hpp

#include <Fcn1D/Fcn1D.hpp>
#include <bininfo.hpp>
#include <TRandom3.h>


namespace peak_search
{

/// @brief Given a 'histo_1D_t' object, a PDF, and a TOTAL number of events to generate 'N', return toy data for each bin (poisson-generated)
/// @param histogram histogram object to put toy events in. all current bin contents will be deleted! 
/// @param PDF PDF to generate with. (average) total events to be generated is integral of PDF from xmin to xmax times N. 
/// @param N let X := int_x0^x1 PDF(x)dx. In this histogram, the average total number of events generated will be <N_hist> = X*N. 
/// @param generator random number generator to use. 
void generate_toy_events(histo_1D_t& histogram, Fcn1D* PDF, double N, TRandom3& generator);

};


#endif