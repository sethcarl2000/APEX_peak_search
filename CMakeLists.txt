#---------------------------------------------------------
#   I copied most of this CMakeLists.txt from: https://cliutils.gitlab.io/modern-cmake/examples/root-dict/README.html 
cmake_minimum_required(VERSION 3.11)

#---------------------------------------------------------
#   Define the project 
project(ApexPeakSearchLibrary VERSION 0.1
        LANGUAGES CXX)  

#---------------------------------------------------------
#   Define this preprocessor macro to enable debug statements 
option(DEBUG "Enable debug printing" OFF)
if (DEBUG)
    add_compile_definitions(DEBUG)
endif()

#---------------------------------------------------------
#   Which include directories to use. 
set(CMAKE_CXX_STANDARD
    20
    CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

#---------------------------------------------------------
#   Find ROOT (necessary)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO Physics MathCore)
#   Find Eigen3 
find_package(Eigen3 5.0 REQUIRED)

#---------------------------------------------------------
#   Define sources and headers
set(sources 
    ExponentialPoly.cpp
    log_likelihood.cpp
    gauss_integrate.cpp
    fit_exponential_poly.cpp
    numbers.cpp
    bininfo.cpp
    best_likelihood_fit.cpp
    chisquare.cpp
)

set(headers 
    ExponentialPoly.hpp
    log_likelihood.hpp
    gauss_integrate.hpp
    fit_exponential_poly.hpp
    FitResult.hpp
    pdf_fcn.hpp
    numbers.hpp 
    bininfo.hpp
    best_likelihood_fit.hpp
    chisquare.hpp
)

#---------------------------------------------------------
#   Generate the ROOT dictionary
root_generate_dictionary(G__ApexPeakSearch ${headers} LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
add_library(ApexPeakSearch SHARED ${sources} ${headers} G__ApexPeakSearch.cxx)
target_include_directories(ApexPeakSearch PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(ApexPeakSearch PUBLIC ROOT::Core Eigen3::Eigen)

