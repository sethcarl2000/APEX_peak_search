#include <GP.hpp>
#include <Fcn1D/Gauss.hpp>
// ROOT headers
#include <TGraph.h> 
#include <TGraphErrors.h>
#include <TCanvas.h> 
#include <TRandom3.h> 
#include <TError.h> 
// stdlib headers
#include <vector>

void test_draw_GP(int n_pts=50)
{
    using namespace peak_search; 

    std::vector<GP::Point> points; points.reserve(n_pts+2); 

    // Error 
    const double standard_variance = 0.05; 

    TRandom3 mrand; 
    points.emplace_back(-0.050, 0.0 + mrand.Gaus()*std::sqrt(standard_variance), standard_variance); 
    points.emplace_back(-0.025, 0.025 + mrand.Gaus()*std::sqrt(standard_variance), standard_variance);   

    Gauss kernel(2., 0., 1.5); 

    std::vector<double> 
        pts_x(n_pts,0), 
        pts_actual(n_pts,0.), pts_actual_error(n_pts,0.); 

    double x =0.;
    double dx =1./(10.); 
    for (int i=0; i<n_pts; i++) {
        

        if (i == 25) {
            x += dx*14.; 
        }

        std::vector<GP::Point> new_pt_v{{x, 0., 0.}}; 
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-0ppppppp
        
        auto& new_pt = new_pt_v[0]; 

        

        GP::Compute(points, new_pt_v, kernel); 


        std::printf("new pt: (x, y, variance): (%+3.1f, %+3.1f, %+3.1f)\n", new_pt.x, new_pt.y, new_pt.variance); 

        if (double yy = new_pt.y; yy != yy) {
            Error(__func__, "NaN encountered; pt. %i/%i, x=%.3f", i,n_pts-1, new_pt.x); 
            n_pts = i; 
            break; 
        }
        
        // pick new point
        pts_x[i] = new_pt.x; 

        new_pt.y += std::sqrt(new_pt.variance + standard_variance)*mrand.Gaus(); 
        new_pt.variance = standard_variance; 

        pts_actual[i] = new_pt.y; 
        pts_actual_error[i] = std::sqrt(standard_variance); 
        points.emplace_back(new_pt);

        x += dx; 
    }

    //now, use these points to 'predict' the mean 
    std::vector<GP::Point> points_predict;

    x = points.front().x; 
    do {
        points_predict.emplace_back(x, 0., standard_variance); 
        x += dx; 
    } while (x < points.back().x-dx/2.); 

    GP::Compute(points, points_predict, kernel); 

    std::vector<double> 
        pts_x_predict,
        pts_predict,
        pts_predict_error; 

    pts_x_predict.reserve(points_predict.size()); 
    pts_predict.reserve(points_predict.size());
    pts_actual_error.reserve(points_predict.size()); 

    for (const auto& pt : points_predict) { 
        pts_x_predict.emplace_back(pt.x); 
        pts_predict.push_back(pt.y); 
        pts_predict_error.push_back(std::sqrt(pt.variance + standard_variance)); 
    }

    new TCanvas; 
    auto g_pred_error = new TGraphErrors(pts_x_predict.size(), pts_x_predict.data(), pts_predict.data(), nullptr, pts_predict_error.data()); 
    g_pred_error->SetFillColor(kGray); 
    g_pred_error->SetLineStyle(0); 
    g_pred_error->Draw("A 3"); 

    auto g_predict = new TGraph(pts_x_predict.size(), pts_x_predict.data(), pts_predict.data()); 
    g_predict->SetLineStyle(kSolid); 
    g_predict->SetLineColor(kBlack);
    g_predict->SetLineWidth(2);  
    g_predict->Draw("SAME L"); 

    auto g = new TGraphErrors(n_pts, pts_x.data(), pts_actual.data(), nullptr, pts_actual_error.data()); 
    g->SetMarkerStyle(kPlus); 
    g->Draw("P Z"); 

    
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++-----------------------------------b ); 

    return; 
}