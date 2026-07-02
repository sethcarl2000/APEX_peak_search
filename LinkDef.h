#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// add each user defined class here, with the syntax: 
//#pragma link C++ class MyClass+; 

#pragma link C++ struct peak_search::Fcn1D+; 
#pragma link C++ struct peak_search::ExponentialPoly+; 
#pragma link C++ struct peak_search::LegendrePoly+; 
#pragma link C++ struct peak_search::Gauss+; 
#pragma link C++ struct peak_search::FcnSum+; 

#endif 