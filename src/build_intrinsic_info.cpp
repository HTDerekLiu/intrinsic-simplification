#include "build_intrinsic_info.h"

void build_intrinsic_info(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs)
{
    build_glue_map(F,G);
    build_edge_lengths(V,F,l);
    build_angular_coordinates(F,G,l,A,v2fs);
}