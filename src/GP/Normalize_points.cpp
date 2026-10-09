#include <GP.hpp>
// stdlib
#include <cmath> 

namespace peak_search
{
namespace GP
{

std::vector<Point> Normalize_points(const std::vector<Point>& data)
{
    if (data.empty()) return {}; 

    std::vector<GP::Point> output; output.reserve(data.size()); 

    // for the given field, scale all values so they fall in the range [-1, +1]
    auto get_maxmin = [&data](double Point::*field, double& min, double& max) {
        min=+1e30; max=-1e30; 
        for (const auto& pt : data) { 
            min = std::min(pt.*field, min); 
            max = std::max(pt.*field, max); 
        }
    };

    double xmin, xmax; 
    get_maxmin(&Point::x, xmin, xmax); 

    double ymin, ymax; 
    get_maxmin(&Point::y, ymin, ymax); 

    double y_mean  = (ymax + ymin)/2.;
    double y_scale = (ymax - ymin)/2.; 

    double x_mean  = (xmax + xmin)/2.;
    double x_scale = (xmax - xmin)/2.;

    std::vector<Point> points; points.reserve(data.size()); 
    
    for (const auto& pt : data) {

        double x = pt.x; 
        //scale the bin center 
        x = (x - x_mean)/x_scale; 

        // find the log of the bin contents (and the leading central moment of the uncertainty) 
        // (we're going to fit in log-space)
        double y = pt.y; 
        double variance = pt.variance / (y_scale*y_scale); 

        y = (y - y_mean)/y_scale; 
    
        points.emplace_back(x, y, variance);
    }
    return points; 
}

}
}