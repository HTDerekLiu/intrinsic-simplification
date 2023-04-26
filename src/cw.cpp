#include "cw.h"
Eigen::Vector2i cw(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs)
{
    return next(twin(G, fs));
}