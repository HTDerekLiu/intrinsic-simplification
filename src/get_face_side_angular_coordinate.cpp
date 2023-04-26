#include "get_face_side_angular_coordinate.h"

double get_face_side_angular_coordinate(
    const Eigen::MatrixXd & A,
    const Eigen::Vector2i & fs)
{
    using namespace Eigen;
    using namespace global_variables;

    if (is_same_face_side(fs, GHOST_FACE_SIDE))
        return DOUBLE_INF;
    else
        return A(fs(0), fs(1));
}