#include "save_matrix.h"

// Stolen from https://github.com/nmwsharp/nonmanifold-Laplacian src/main.cpp
void save_matrix(
                 Eigen::SparseMatrix<double>& matrix,
                 std::ostream& out)
{
    out << std::setprecision(16);

    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it ) {
            double val       = it.value();
            size_t iRow = it.row();
            size_t iCol = it.col();

            out << (iRow + 1) << " " << (iCol + 1) << " " << val << std::endl;
        }
    }
}

void save_matrix(
                Eigen::SparseMatrix<double>& matrix,
                std::string filename)
{
    // WARNING: this follows matlab convention and thus is 1-indexed

    std::cout << "Writing sparse matrix to: " << filename << std::endl;

    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("failed to open output file " + filename);
    }
    save_matrix(matrix, outFile);
    outFile.close();
}
