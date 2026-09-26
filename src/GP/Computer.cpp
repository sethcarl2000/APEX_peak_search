#include <GP.hpp>

namespace peak_search
{
namespace GP
{

void Computer::Compute(std::vector<Point> inputs, std::vector<Point>& outputs) const
{
    //scale points
    double var_scale = 1./(y_scale*y_scale); 

    for (auto& pt : inputs) {
        pt.x = (pt.x - x_mean)/x_scale;
        pt.y = (pt.y - y_mean)/y_scale; 

        pt.variance *= var_scale; 
    }

    for (auto& pt : outputs) {
        pt.x = (pt.x - x_mean)/x_scale;
        pt.y = (pt.y - y_mean)/y_scale; 

        pt.variance *= var_scale; 
    }

    GP::Compute(inputs, outputs, fKernel); 

    //now, de-scale the y-points
    for (auto& pt : outputs) {
        pt.x = (x_scale*pt.x) + x_mean;
        pt.y = (y_scale*pt.y) + y_mean;

        pt.variance /= var_scale; 
    }
}

}
}