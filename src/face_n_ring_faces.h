#ifndef FACE_N_RING_FACES
#define FACE_N_RING_FACES

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <iostream>

#include "twin.h"
#include "next.h"
#include "is_same_face_side.h"
#include "global_variables.h"

/*
    this function returns the neighboring n-ring faces of a given face

    Inputs
    G: |F|x3x2 array of gluing map 
    f_src: index of the center face
    num_rings: int of the number of rings to extract

    Outputs
    fIdx: list of neighboring face indices
*/
void face_n_ring_faces(
    const Eigen::MatrixXi & G,
    const int & f_src,
    const int & num_rings,
    Eigen::VectorXi & fIdx);
#endif