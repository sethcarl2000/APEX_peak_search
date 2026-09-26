#include <GP.hpp>
// ROOT
#include <TGraph.h>
#include <TStyle.h> 
#include <TPad.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TError.h>
#include <TH1F.h> 
// Eigen
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
// stdlib
#include <vector> 
#include <cstdio> 
#include <iostream> 


using Eigen::VectorXd, Eigen::MatrixXd, Eigen::SelfAdjointEigenSolver; 

using namespace peak_search; 

/// @brief Given x points and kernel, generate random Y-values from given kernel (with zero mean)
/// @param X input x-values
/// @param kernel kernel to use 
/// @return Y-values randomly generated
VectorXd Gen_GP(VectorXd& X, const GP::Kernel& kernel); 

//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
void draw_random_GP() 
{ 
    std::vector<double> pars{ 1., 5. }; 
    
    GP::Kernel kernel{
        [](double x1, double x2, const std::vector<double>& params){
            
            if (params.size() != 2) return peak_search::numbers::nan; 
            double amplitude{params[0]}, length_scale{params[1]}; 
            double arg = (x1 - x2)/length_scale; 
            return amplitude * std::exp( -0.5*arg*arg ); 
        }, pars
    };

    if (!gPad) { 
        new TCanvas; 
    } else {
        gPad->Clear(); 
    }

    const int n_pts = 20; 

    VectorXd X(n_pts); 

    for (int i=0; i<n_pts; i++) X(i) = (double)i; 

    VectorXd Y = Gen_GP(X, kernel);

    for (int i=0; i<n_pts; i++) {
        //std::printf(" %.4f, %+.4f\n", X(i), Y(i));
    }
    auto frame = gPad->DrawFrame(-0.5,-2.5, ((double)n_pts)-0.5,+2.5); 
    frame->SetTitle(";x_{i};y_{i}"); 
    auto graph = new TGraph(n_pts, X.data(), Y.data()); 
    graph->SetMarkerStyle(kOpenCircle); 
    graph->Draw("SAME LP"); 

    gPad->Modified(); 
    gPad->Update(); 

    //canv->SaveAs("plots/gp_example_2.gif+100"); 
    //canv->Clear(); 

    return; 
}
//________________________________________________________________________________________________________
VectorXd Gen_GP(VectorXd& X, const GP::Kernel& kernel)
{
    static TRandom3 myrand; 

    const size_t n_pts = X.size(); 

    MatrixXd cov(n_pts, n_pts);

    for (int i=0; i<n_pts; i++) 
        for (int j=i; j<n_pts; j++) { cov(i,j) = cov(j,i) = kernel(X(i), X(j)); }

    //now, solve the cov. matrix. 
    SelfAdjointEigenSolver<MatrixXd> es(cov);

    std::cout << cov << "\n"; 

    VectorXd Y = VectorXd::Zero(n_pts); 
    
    if (es.info() != Eigen::Success) {
        Error(__func__, "Something went wrong fiding the eigenvectors of the cov. matrix"); 
        return Y; 
    }

    MatrixXd eVecs = es.eigenvectors(); 
    VectorXd eVals = es.eigenvalues(); 
    int n_frames = 20; 

    VectorXd E(n_pts);
    for (int i=0; i<n_pts; i++) {

        double gaus = myrand.Gaus(); 

        for (int j=0; j<n_pts; j++) {
            E(j) = gaus*std::sqrt(std::fabs(eVals(i)))*eVecs(i,j); 
        }
        Y += E; 
    }
    return Y; 
} 
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________