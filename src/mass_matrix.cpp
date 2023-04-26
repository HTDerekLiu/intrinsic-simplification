#include "mass_matrix.h"

void mass_matrix(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::SparseMatrix<double> & M)
{
    using namespace Eigen;

    VectorXd VA;
    vertex_areas(F, l, VA);

    int nV = F.maxCoeff()+1;
    M.resize(nV,nV);
    M = VA.asDiagonal();
}