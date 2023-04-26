#include "opposite_corner_angle.h"

double opposite_corner_angle(
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs)
{
    using namespace global_variables;
    using namespace std;
    
    // extract edge lengths
    double lij, ljk, lki;
    lij = l(fs(0), fs(1));
    ljk = l(fs(0), (fs(1)+1) % 3);
    lki = l(fs(0), (fs(1)+2) % 3);

    // uniformly scale the edge to increase numerics
    double inv_lsum = 1.0 / (lij + ljk + lki);
    lij *= inv_lsum; // equiv to lij /= lsum
    ljk *= inv_lsum; // equiv to ljk /= lsum
    lki *= inv_lsum; // equiv to lki /= lsum

    // law of cosines
    double d = (ljk*ljk + lki*lki - lij*lij) / (2*ljk*lki);

    // compute arccos
    if (abs(d) <= 1) 
    {
        // satisfy triangle inequality
        return acos(d);
    }
    else if (abs(d)>1 && abs(d)<(1+EPS))
    {
        // slightly violate triangle inequality, clip it
        d = std::max(-1.0, std::min(d, 1.0)); // clip to -1~1
        return acos(d);
    }
    else
    {
        std::cout << "WARNING: In opposite_corner_angle.cpp, the triangle violates triangle inequality" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
}
