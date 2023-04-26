#include "twin.h"

Eigen::Vector2i twin(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs)
{
    return Eigen::Vector2i(G(fs(0), fs(1)*2), G(fs(0), fs(1)*2+1));
}
