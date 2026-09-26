#include <GP.hpp>
// stdlib
#include <cmath> 

namespace peak_search
{
namespace GP
{

std::vector<Point> Normalize_points(const Histo1D& data)
{
    // for the given field, scale all values so they fall in the range [-1, +1]
    auto get_maxmin = [&data](double HistoBin::*field, double& min, double& max) {
        min=+1e30; max=-1e30; 
        for (const auto& bin : data.bins) { 
            min = std::min(bin.*field, min); 
            max = std::max(bin.*field, max); 
        }
    };

    double xmin{data.GetXmin()}, xmax{data.GetXmax()}; 

    double ymin, ymax; 
    get_maxmin(&HistoBin::N, ymin, ymax); 

    ymin = std::log(ymin); 
    ymax = std::log(ymax); 

    double y_mean  = (ymax + ymin)/2.;
    double y_scale = (ymax - ymin)/2.; 

    double x_mean  = (xmax + xmin)/2.;
    double x_scale = (xmax - xmin)/2.;

    std::vector<Point> points; points.reserve(data.GetNbins()); 
    
    for (const auto& bin : data.bins) {

        //find the bin center
        double x = (bin.xmin + bin.xmax)/2.;
        
        //scale the bin center 
        x = (x - x_mean)/x_scale; 

        // find the log of the bin contents (and the leading central moment of the uncertainty) 
        // (we're going to fit in log-space)
        double y = std::log(bin.N); 
        double stddev = 1./std::sqrt(bin.N);
        
        y = (y - y_mean)/y_scale; 
        stddev *= 1./y_scale; 
    
        points.emplace_back(x, y, stddev*stddev);

    }
    return points; 
}

}
}