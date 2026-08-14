
#include <read_model_from_file.hpp>

//stdlib headers
#include <fstream> 
#include <stdexcept> 
#include <sstream> 
#include <fstream> 
#include <vector> 
//ROOT headers
#include <TString.h> 
#include <TError.h> 



namespace peak_search
{

void read_model_from_file(const std::string& file_path, Fcn1D* model) 
{
    //check if the model passed is null
    if (!model) {
        throw std::invalid_argument(Form(
            "in <%s>: model ptr is null, must be initialized before being passed.",
            __func__
        ));
        return; 
    }

    //try to open the file
    std::ifstream infile(file_path); 

    if (!infile.is_open()) {
        throw std::runtime_error(Form(
            "in <%s>: Failed to open input file: '%s'.",
            __func__, file_path.c_str()
        ));
        return;
    }

    //get & reset the list of coefficients
    auto& coefficients = model->GetParams(); 
    coefficients.clear(); 

    std::string line; 
    int elems_read=0;

    int line_num=0; 
    //read each line of the file
    while (std::getline(infile, line)) {

        ++line_num; 

        //break the line up into tokens
        std::istringstream iss(line); 

        std::string token; 

        //check the first token of the line. if it's an '#' symbol, then 
        iss >> token; 
        if (token.empty() || (token[0]=='#')) continue; 

        //now, try converting this token to an int
        int elem_num=0; 
        double elem_val=0; 

        try {

            elem_num = std::stoi(token); 
        
        } catch (const std::exception& e) {

            Warning(__func__, "<file '%s', line %i>: exception caught trying to convert string '%s' to integer. Unable to parse line.\n what(): %s\nLine skipped.", 
                file_path.c_str(), line_num,
                token.c_str(), e.what()
            );
            continue; 
        }

        //check to make sure the element index is right. 
        if (elems_read != elem_num) {
            Error(__func__, "<file '%s', line %i>: Expected element index '%i' to be next, but this element has index '%i' (out-of-order?). "
                "File parsing terminated.", 
                file_path.c_str(), line_num, 
                elems_read, elem_num
            ); 
            break; 
        }
        
        iss >> token; 

        try {

            elem_val = std::stod(token); 
        
        } catch (const std::exception& e) {

            Warning(__func__, "<file '%s', line %i>: exception caught trying to convert string '%s' to double. Unable to parse line.\n what(): %s\nLine skipped.", 
                file_path.c_str(), line_num,
                token.c_str(), e.what()
            );
            continue; 
        }

        coefficients.push_back(elem_val); 
        ++elems_read; 
    }   

    infile.close(); 
}

};