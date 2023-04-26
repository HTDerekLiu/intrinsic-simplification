#include "face_areas.h"
void face_areas(
    const Eigen::MatrixXd & l,
    Eigen::VectorXd & a)
{
    using namespace Eigen;

    int nF = l.rows();
    a.resize(nF);

    // Heron's rule
    a = (l.col(0) + l.col(1) + l.col(2)) / 2;
    a = a.array() * (a-l.col(0)).array() * (a-l.col(1)).array() * (a-l.col(2)).array();
    a = a.cwiseSqrt();
}
