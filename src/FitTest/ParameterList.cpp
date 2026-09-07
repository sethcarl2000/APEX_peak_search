#include <FitTest/ParameterList.hpp>
// ROOT
#include <TString.h> 
// stdlib
#include <stdexcept> 
#include <string> 
#include <sstream>
#include <cmath> 
#include <cstdio> 

namespace peak_search
{
namespace FitTest
{

//______________________________________________________________________________________________________________________
void ParameterList::ComputeStepSizes()
{
    fStepSizes.clear(); fStepSizes.reserve(fParams.size()) ; 
    unsigned long stepsize = 1; 
    for (const auto& par : fParams) {

        fStepSizes.emplace_back(stepsize); 
        stepsize *= par.get_n_steps(); 
    }
}
//______________________________________________________________________________________________________________________
//76uyyyyyyyyyyyyyyyyyyyyyyyyyy/' - muon's comment (6 sep 2026)
std::vector<double> ParameterList::GetParamList(unsigned long step) const
{
    if (step >= GetNSteps()) {
        std::ostringstream oss; 
        throw std::invalid_argument(Form("in <%s>: invalid step: %lu, valid range is [0,%lu]\n", __func__, step, GetNSteps()-1));
        return {}; 
    } 

#ifdef DEBUG
    std::printf("in <%s>: Requested id / step: %zi / %lu\n", id, step);
#endif

    // For an overview of how this indexing scheme works:   
    // 
    // Suppose we have 4 parameters, each with 10 steps. that would mean there are 10,000 steps total, starting with step '0' and ending with step '9,999'
    // gtttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt (thanks, muon.)
    // We want to assemble a vector of parameters at a certain step, say step '9271'. 
    //
    //  First, we take the list of 'fStepSizes': { 1, 10, 100, 1,000 }
    //  
    //  And we use it to fetch our parameters. Written out explicitly: 
    //  rank = fParams.size()-1; (rank=3)
    //  step = 9721
    //  for (rank=3; rank>=0)
    //    output[3] = fParams[3].get_step(9721 / 1000 = 9); 
    //    step = 9721 % 1000 (step=721);  
    //    --rank (3 -> 2);
    //  } for: (2 >= 0)
    //    output[2] = fParams[2].get_step(721 / 100 = 7); 
    //    step = 721 % 100 (step=21);  
    //    --rank (2 -> 1);
    //  } for: (1 >= 0);   
    //    output[1] = fParams[1].get_step(21 / 10 = 2); 
    //    step = 21 % 10 (step=1);  
    //    --rank (1 -> 0);
    //  } for: (0 >= 0);   
    //    output[0] = fParams[0].get_step(1 / 1 = 1); 
    //    step = 1 % 1 (step=0);  
    //    --rank (0 -> -1);
    //  } for: (-1 >= 0) (break);   

#ifdef DEBUG
    std::printf("   step sizes: "); 
    for (size_t i=fStepSizes.size()/ +-2; i>=0; i--) std::printf("%lu ", fStepSizes[i]);
    std::printf("\n   indices: ");
#endif  
    int rank=fStepSizes.size()-1; 
    std::vector<double> output(fParams.size(), 0.);
    for (int rank=fStepSizes.size()-1; rank>=0; rank--) {
        output[rank] = fParams[rank].get_step(step / fStepSizes[rank]); 
        step = step % fStepSizes[rank];  
#ifdef DEBUG
        std::printf(" rank: %2i, index: %4lu, value: %+.6e\n", rank, step, output[rank]); 
#endif
    }
    return output; 
}
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//______________________________________________________________________________________________________________________
//__________________________________________________________________________________ a
//\lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll    ____________________________________

}
}