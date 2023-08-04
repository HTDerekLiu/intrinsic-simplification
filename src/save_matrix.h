#ifndef SAVE_MATRIX
#define SAVE_MATRIX

#include <Eigen/Sparse>

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream> // std::ofstream

/*
This function saves a sparse matrix to a file.

Sparse matrices are is exported as a list of triplets. Explicitly, each line of the output file contains the row index, column index, and value of some entry of the matrix. These indices are 1-indexed to make it easy to load in [Matlab](https://www.mathworks.com/help/matlab/ref/spconvert.html).

Inputs:
    matrix: the matrix to output
    filename: a file to save the matrix to (e.g. matrix.spmat)
      or
    out: an std::ostream to write the matrix to
*/
void save_matrix(
                 Eigen::SparseMatrix<double> & matrix,
                 std::string filename);

void save_matrix(
                 Eigen::SparseMatrix<double> & matrix,
                 std::ostream & out);

#endif
