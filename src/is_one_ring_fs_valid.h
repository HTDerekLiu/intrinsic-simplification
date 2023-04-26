#ifndef IS_ONE_RING_FS_VALID
#define IS_ONE_RING_FS_VALID

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "global_variables.h"
#include "get_face_side_vertices.h"
#include "vertex_one_ring_face_sides.h"

bool is_one_ring_fs_valid(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v);

#endif
