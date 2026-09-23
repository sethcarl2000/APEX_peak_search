#include <GP/Point.hpp>
#include <GP/Compute.hpp>
// nlopt
#include <nlopt.h> 
// Eigen
#include <eigen3/Eigen/Core> 
#include <eigen3/Eigen/Dense>
// stdlib
#include <sstream>
#include <string>
#include <stdexcept> 
#include <vector> 
#include <iostream> 

namespace peak_search
{
namespace GP
{

void Compute(const std::vector<Point>& inputs, std::vector<Point>& outputs, const Fcn1D& kernel)
{
    auto wrapper = [&kernel](double x1, double x2){ return kernel(x1-x2); };
    Compute(inputs, outputs, wrapper); 
}
//_____________________________________________________________________________________________________
void Compute(const std::vector<Point>& inputs, std::vector<Point>& outputs, const std::function<double(double,double)>& kernel)
{
    if (outputs.empty()) {
        std::ostringstream oss; 
        oss << "in <"<<__func__<<">: No output points given."; 
        throw std::invalid_argument(oss.str()); 
        return; 
    }

    const int n_inputs = inputs.size(); 
    if (n_inputs < 2) {
        std::ostringstream oss; 
        oss << "in <"<<__func__<<">: Number of input points ("<<n_inputs<<") invalid; must be at least 2."; 
        throw std::invalid_argument(oss.str()); 
        return; 
    }

#ifdef DEBUG
    std::cout << "input points: ~~~~~~~~~~~~~\n"; 
    for (const auto& pt : inputs) {
        std::cout << pt.x << "  " << pt.y << "  " << pt.variance << "\n"; 
    }
#endif

    // construct the B-matrix 
    using Eigen::MatrixXd, Eigen::VectorXd, Eigen::LLT;

    // compute the upper-triangle of the matrix (including diagonal)
    MatrixXd B(n_inputs, n_inputs); 
    VectorXd Y(n_inputs); 

    for (int i=0; i<n_inputs; i++) {

        double xi = inputs[i].x;   
        Y(i)      = inputs[i].y; 

        B(i,i) = kernel(xi,xi) + inputs[i].variance; 
            
        for (int j=i+1; j<n_inputs; j++) {

            double xj = inputs[j].x; 
            B(i,j) = kernel(xi,xj); 
#ifdef DEBUG
            std::cout << "k("<<xi<<"-"<<xj<<") = " << B(i,j) << "\n";
#endif
        }
    }
    // now, symmetrize the matrix
    for (int i=1; i<n_inputs; i++) 
        for (int j=0; j<i; j++) 
            B(i,j) = B(j,i); 

    //now, solve the system 
    Eigen::LDLT L = B.ldlt(); 

#ifdef DEBUG
    std::cout << "matrix: ~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << B << "\n"; 
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"; //*/ 
#endif

    for (auto& pt : outputs) {

        double x = pt.x; 

        VectorXd K(n_inputs); 
        for (int j=0; j<n_inputs; j++) K(j) = kernel(x,inputs[j].x);

        VectorXd J = L.solve(K); 

        pt.y = J.dot(Y); 

        pt.variance = kernel(x,x) - J.dot(K); 
    }
}

}
}