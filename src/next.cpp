#include "next.h"

Eigen::Vector2i next(
    const Eigen::Vector2i & fs)
{
    using namespace global_variables;

    int f = fs(0);
    int s = fs(1);
    if (f != GHOST_INDEX && s != GHOST_INDEX) // if not ghost fs
        return Eigen::Vector2i(f, fast_mod(s+1, 3));
    else
        return fs;
}
