#ifndef IS_GLUE_MAP_VALID
#define IS_GLUE_MAP_VALID

#include <Eigen/Core>
#include <Eigen/Dense>
#include "global_variables.h"
#include "build_glue_map.h"

bool is_glue_map_valid(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G);

#endif