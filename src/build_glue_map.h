#ifndef BUILD_GLUE_MAP
#define BUILD_GLUE_MAP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>

#include <igl/sortrows.h>

#include "get_face_side_vertices.h"
#include "global_variables.h"
#include "glue_face_sides.h"

/*
Build the glue map for a triangle mesh possibly with boundary.

Inputs
    F: |F|x3 matrix of face list
Outputs
    G: |F|x6 array of gluing map. 

Notes    
For each face side (f, s):
- G[f,2*s] stores which face f' it glues to
- G[f,2*s+1] stores which side s' it glues to 
Thus (f,s) and (f',s') are glued together. If f' or s' is GHOST_INDEX (defined in global_variables.h) it means that (f,s) is a boundary face side and (f',s') does not exist. 
*/

void build_glue_map(
    const Eigen::MatrixXi & F,
    Eigen::MatrixXi & G);
#endif