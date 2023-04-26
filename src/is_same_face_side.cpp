#include "is_same_face_side.h"
bool is_same_face_side(
    const Eigen::Vector2i & fs0,
    const Eigen::Vector2i & fs1)
{
    return (fs0(0)==fs1(0) && fs0(1)==fs1(1));
}