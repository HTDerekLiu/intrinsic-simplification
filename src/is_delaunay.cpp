#include "is_delaunay.h"
bool is_delaunay(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs)
{
    using namespace global_variables;

    if (is_boundary_face_side(G,fs)) // this is boundary fs
    {
        double theta = opposite_corner_angle(l, fs);
        return (theta <= (M_PI + EPS));
    }
    else // interoir fs
    {
        double theta_A = opposite_corner_angle(l, fs);
        double theta_B = opposite_corner_angle(l, twin(G,fs));
        return ((theta_A + theta_B) <= (M_PI + EPS));
    }

}