#ifndef MG_PRECOMPUTE
#define MG_PRECOMPUTE

#include <vector>
#include <iostream>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <mg_data.h>
#include <get_prolongation.h>
#include <coarsen_mesh.h>

// Build up a multi-level hierarchy for the surface multigrid solver
//
// Inputs:
//   ratio  decimation rate between levels
//   nVCoarset  desired number of output faces for the coarest level
//   dec_type decimation type (0:qslim, 1:midpoint, 2:vertex removal)
//
// Outputs:
//   mg    a vector that contains precomputed information for each level

void mg_precompute(
    const Eigen::MatrixXi & Fi,
    const Eigen::MatrixXi & Gi,
    const Eigen::MatrixXd & li,
    const Eigen::MatrixXd & Ai,
    const Eigen::MatrixXi & v2fsi, 
    const double & ratio,
    const int & nV_coarsest,
    const double & mixture_weight,
    std::vector<mg_data> & mg);

void mg_precompute(
    const Eigen::MatrixXi & Fi,
    const Eigen::MatrixXi & Gi,
    const Eigen::MatrixXd & li,
    const Eigen::MatrixXd & Ai,
    const Eigen::MatrixXi & v2fsi, 
    const double & ratio,
    const int & nV_coarsest,
    std::vector<mg_data> & mg);
#endif