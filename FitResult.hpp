#ifndef peak_search_FitResult_hpp
#define peak_search_FitResult_hpp

#include <TError.h> 

namespace peak_search
{

namespace Status {
    enum bit : int { 
        kNull    =0, 
        kFail    =-1, 
        kSuccess =+1
    };
}

//generic fit-result template, to allow passing of fit-result data 
template<typename FitDataType> struct FitResult {

    //check to make sure the fit-result data type is default-constructible
    static_assert(std::is_default_constructible_v<FitDataType>, "FitDataType template param is not default-constructible");

    
    FitDataType data; 
    Status::bit status{Status::kNull};

    inline explicit operator bool() const { return status == Status::kSuccess; }

    inline explicit operator FitDataType() const { 
        return data; 
    }

    //null result
    static FitResult<FitDataType> Null() { return FitResult<FitDataType>{ .status = Status::kNull; }; }
    //failed result
    static FitResult<FitDataType> Fail() { return FitResult<FitDataType>{ .status = Status::kFail; }; }
};


};
#endif