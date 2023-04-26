#include "ccw.h"

Eigen::Vector2i ccw(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs)
{
    return twin(G, next(next(fs)));
}