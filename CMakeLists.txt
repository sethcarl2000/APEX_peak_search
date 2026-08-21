#---------------------------------------------------------
#   I copied most of this CMakeLists.txt from: https://cliutils.gitlab.io/modern-cmake/examples/root-dict/README.html 
cmake_minimum_required(VERSION 3.16)

#---------------------------------------------------------
#   Define the project 
project(ApexPeakSearchLibrary VERSION 0.1
        LANGUAGES CXX)  

#---------------------------------------------------------
#   Find ROOT (necessary)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO Physics MathCore)
#   Find Eigen3 
find_package(Eigen3 5.0 REQUIRED)

#---------------------------------------------------------
#   Define sources and headers
set(sources 
    src/legendre_polynomial.cpp
    src/log_likelihood.cpp
    src/gauss_integrate.cpp
    src/fit_exponential_poly.cpp
    src/fit_exponential_legendre.cpp
    src/numbers.cpp
    src/best_likelihood_fit.cpp
    src/chisquare.cpp
    src/compute_Q0.cpp
    src/newton_optimizer.cpp
    src/read_model_from_file.cpp
    src/generate_toy_events.cpp
    src/make_histogram_copy.cpp
    src/FitTestKernel.cpp
    src/FitTestThreadManager.cpp
    # Fcn1D 
    src/Fcn1D/Fcn1D.cpp
    src/Fcn1D/ExponentialPoly.cpp
    src/Fcn1D/ExponentialLegendre.cpp
    src/Fcn1D/LegendrePoly.cpp
    src/Fcn1D/Gauss.cpp
    src/Fcn1D/FcnSum.cpp
)

set(headers 
    include/legendre_polynomial.hpp
    include/log_likelihood.hpp
    include/gauss_integrate.hpp
    include/fit_exponential_poly.hpp
    include/fit_exponential_legendre.hpp
    include/FitResult.hpp
    include/numbers.hpp 
    include/best_likelihood_fit.hpp
    include/chisquare.hpp
    include/fit_parameter.hpp
    include/compute_Q0.hpp
    include/newton_optimizer.hpp
    include/read_model_from_file.hpp
    include/generate_toy_events.hpp
    include/make_histogram_copy.hpp
    include/HistoBin.hpp
    include/Histo1D.hpp
    include/FitTestFunction.hpp
    include/FitTestKernel.hpp
    include/FitTestThreadManager.hpp
    include/SignalFit.hpp
    # Fcn1D 
    include/Fcn1D/Fcn1D.hpp
    include/Fcn1D/ExponentialPoly.hpp
    include/Fcn1D/ExponentialLegendre.hpp
    include/Fcn1D/LegendrePoly.hpp
    include/Fcn1D/Gauss.hpp
    include/Fcn1D/FcnSum.hpp
)

#---------------------------------------------------------
#   Generate the ROOT dictionary
#root_generate_dictionary(G__ApexPeakSearch ${headers} LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries

#
#---------------------------------------------------------
#   Define this preprocessor macro to enable debug statements 


add_library(ApexPeakSearch SHARED ${sources} )# G__ApexPeakSearch.cxx)

#-----------------------------------------------------------------------
#   Add some compile options
option(DEBUG "Enable debug printing in all code" OFF)
if (DEBUG)
    target_compile_definitions(ApexPeakSearch PRIVATE DEBUG)
endif()
#
option(DEBUG_NEWTON "Debug printing for 'newton_optimizer'" OFF)
if (DEBUG_NEWTON)
    target_compile_definitions(ApexPeakSearch PRIVATE DEBUG_NEWTON)
endif()


target_include_directories(ApexPeakSearch PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(ApexPeakSearch PUBLIC ROOT::Core Eigen3::Eigen)

root_generate_dictionary(G__ApexPeakSearch ${headers} LINKDEF LinkDef.h MODULE ApexPeakSearch)