#include "save_complex_matrix.h"

void save_complex_matrix(
                         Eigen::SparseMatrix<std::complex<double>> & matrix,
                         std::string filename) {
    std::string ext = igl::extension(filename);
    if (ext.empty()) {
        save_complex_matrix(matrix, filename + "_re", filename + "_im");
    } else {
        std::string plain_name = filename.substr(0, filename.size() - ext.size() - 1);
        save_complex_matrix(matrix, plain_name + "_re." + ext, plain_name + "_im." + ext);
    }
}

void save_complex_matrix(
                         Eigen::SparseMatrix<std::complex<double>> & matrix,
                         std::string filename_real,
                         std::string filename_imag) {
    // it would be nice to avoid these copies
    Eigen::SparseMatrix<double> real = matrix.real();
    save_matrix(real, filename_real);
    Eigen::SparseMatrix<double> imag = matrix.imag();
    save_matrix(imag, filename_imag);
}
