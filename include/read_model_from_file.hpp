#ifndef peak_search_read_model_from_file_hpp
#define peak_search_read_model_from_file_hpp

#include <Fcn1D/Fcn1D.hpp>
//stdlib headers
#include <string> 

namespace peak_search 
{

/// @brief reads a Fcn1D model from a file. 
/// @param file_path file path, in which a list of coefficients is saved. see models in the 'data/models' directory to see formatting patterns. 
/// @param model model in which to store the coefficients. Pass an already initialzed model! otherwise function fails. 
void read_model_from_file(const std::string& file_path, Fcn1D* model);

};



#endif