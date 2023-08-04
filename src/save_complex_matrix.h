#ifndef SAVE_COMPLEX_MATRIX
#define SAVE_COMPLEX_MATRIX

#include <Eigen/Sparse>

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream> // std::ofstream

#include <igl/extension.h>
#include <save_matrix.h>

/*
  This function saves a complex sparse matrix to a pair of files, encoding its real and imaginary parts as separate, real sparse matrices.

  Sparse matrices are is exported as a list of triplets. Explicitly, each line of the output file contains the row index, column index, and value of some entry of the matrix. These indices are 1-indexed to make it easy to load in [Matlab](https://www.mathworks.com/help/matlab/ref/spconvert.html).

  Inputs:
  matrix: the complex matrix to output
  filename: a base filename to save the matrix to. For example, if filename=matrix.spmat,
            then the real part is saved to matrix_re.spmat, and the imaginary part
            is saved to matrix_im.spmat
    or

  filename_real: the filename for the real part
  filename_imag: the filename for the imaginary part
*/
void save_complex_matrix(
                         Eigen::SparseMatrix<std::complex<double>> & matrix,
                         std::string filename);

void save_complex_matrix(
                         Eigen::SparseMatrix<std::complex<double>> & matrix,
                         std::string filename_real,
                         std::string filename_imag);

#endif
