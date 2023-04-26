#ifndef REMOVE_EAR_VERTEX
#define REMOVE_EAR_VERTEX

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <set>
#include <stdexcept>

#include "is_ear_vertex.h"
#include "flatten_diamond_mesh.h"
#include "compute_barycentric_robust.h"
#include "get_face_side_vertices.h"
#include "twin.h"
#include "next.h"
#include "global_variables.h"
#include "glue_face_sides.h"
#include "flip_to_delaunay_local.h"
#include "cw.h"

/*
    This function removes an ear vertex

    Inputs
    v: the vertex index we try to remove

    V: |V|x3 array of vertex locations
    F: |F|x3 array of face list
    G: |F|x3x2 array of gluing map 
    l: |F|x3 array of edge-lengths array
    A: |F|x3 array of signpost angles
    v2fs: |V|x2 array of face-sides
    BC: |BC|x3 array of barycentric coordinates 
    F2V: |F| length list of lists, where F2V[f] gives you indices in BC.
    T: |V|xch*3 array of tangent vectors

    Outputs
    All the variables are updated in-place

    Notation
                 /  \      
                /    \      
               /      \    
              /        \    
             / fs_twin  \  
    -----v_left ----- v_right -------
             \    fs    /   
              \        /     
               \      /       
                \    /       
                 \  /         
                  v
*/
void remove_ear_vertex(
    const int & v,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V,
    Eigen::MatrixXd & T);

#endif