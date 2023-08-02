#ifndef GET_PROLONGATION
#define GET_PROLONGATION

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>

#include "remove_unreferenced_intrinsic.h"

/*
compute the prolongation operator

Inputs:
    F: (redundant) face list with unreferenced vertices, the output of "coarsen_mesh"
    BC: barycentrinc points (with size equals to V.rows()), also the output of "coarsen_mesh"
    F2V: face to BC indices
    nV: number of fine vertices

Outputs:
    P: |V|x|V_coarse| scalar prolongation matrix

WARNING: 
BC is the output of coarsen_mesh, which has size |nV_fine|x3 (with a lot of redundant rows in it)
*/

void get_prolongation(
	const Eigen::MatrixXi & F, 
    const Eigen::MatrixXd & BC, 
    const std::vector<std::vector<int>> & F2V,
    const std::map<int, int> & IMV,
    Eigen::SparseMatrix<double> & P);

#endif
