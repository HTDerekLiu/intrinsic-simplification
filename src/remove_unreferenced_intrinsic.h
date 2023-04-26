#ifndef REMOVE_UNREFERENCED_INTRINSIC
#define REMOVE_UNREFERENCED_INTRINSIC

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>
#include <vector>
#include <iostream>
#include <igl/slice.h>
#include <map>
#include <unordered_set>
#include "remove_vector_element.h"
#include "global_variables.h"

/*
Remove unreferenced vertex indices, updating F and intrinsic data accordingly

Input:
  F     face list
  optionally other intrinsic data to be updated in place

Outputs:
  RF    reindexed face list
  IMV   index map for vertices such that IMV[v_fine] = v_coarse
  IMF   index map for faces such that IMF[f_fine] = f_coarse
  vIdx  a subset of vertex indices, such that RV = V(vIdx,:) 
  fIdx  a subset of face indices, such that RF = F(fIdx,:)

Note:
  we can get RV by igl::slice(V,vIdx,1,RV);
*/

void remove_unreferenced_intrinsic(
    Eigen::MatrixXi & F,
    std::map<int, int> & IMV,
    std::map<int, int> & IMF,
    Eigen::VectorXi & vIdx,
    Eigen::VectorXi & fIdx);

void remove_unreferenced_intrinsic(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,  
    std::vector<std::vector<int>> & F2V,  
    std::map<int, int> & IMV,
    std::map<int, int> & IMF,
    Eigen::VectorXi & vIdx,
    Eigen::VectorXi & fIdx);

#endif