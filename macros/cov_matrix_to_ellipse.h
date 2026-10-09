#ifndef cov_matrix_to_ellipse_h
#define cov_matrix_to_ellipse_h

// ROOT
#include <TError.h> 
// stdlib
#include <limits> 
#include <cmath> 

/// @brief Given a 2x2 covariance matrix, return the major/minor axes, and angle with x-axis. 
/// @param C11 First diag. element of cov. matrix
/// @param C22 Second diag. element of cov. matrix
/// @param C12 off-diag element of cov. matrix 
/// @param l1 length of first major-axis (horizontal)
/// @param l2 length of second major-axis (vertical)
/// @param theta angle between v1 (first eigenvector) and x-axis
void cov_matrix_to_ellipse(double C11, double C22, double C12, double &l1, double &l2, double &theta)
{
    static const double NaN = std::numeric_limits<double>::quiet_NaN(); 
    l1 = l2 = theta = NaN; 
    
    double trC  = C11 + C22; 
    double detC = C11*C22 - C12*C12; 

    if (trC*trC - 4.*detC < 0.) {
        Error(__func__, 
            "Eigenvalues for following cov matrix are not real:\n"
            "   %+10.4e %+10.4e\n"
            "   %+10.4e %+10.4e\n"
            ,
            C11, C12, 
            C12, C22
        ); 
        return; 
    }

    l1 = 0.5*( trC + std::sqrt(trC*trC - 4.*detC) ); 
    l2 = 0.5*( trC - std::sqrt(trC*trC - 4.*detC) ); 

    if (l1 < 0. || l2 < 0.) {
        Error(__func__, 
            "Eigenvalue(s) for following cov matrix are negative:\n"
            "   %+10.4e %+10.4e\n"
            "   %+10.4e %+10.4e\n"
            " E-vals:\n"
            "   %+10.4e %+10.4e\n"
            ,
            C11, C12, 
            C12, C22,
            l1, l2
        ); 
        l1 = l2 = NaN; 
        return; 
    }

    theta = std::atan( (l1 - C11)/C12 );
    
    return; 
}

#endif