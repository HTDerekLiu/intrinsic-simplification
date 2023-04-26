#ifndef IS_V2FS_VALID
#define IS_V2FS_VALID

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "global_variables.h"
#include "get_face_side_vertices.h"

bool is_v2fs_valid(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & v2fs);

#endif