#include "get_vertex_connection_laplacian.h"

#include "twin.h"
#include "build_angular_coordinates.h"

void get_vertex_connection_laplacian(
                                     const Eigen::MatrixXi & F,
                                     const Eigen::MatrixXi & G,
                                     const Eigen::MatrixXd & l,
                                     Eigen::SparseMatrix<std::complex<double>> & L)
{
    Eigen::MatrixXd theta;
    Eigen::MatrixXi v2fs;
    build_angular_coordinates(F, G, l, theta, v2fs);
    get_vertex_connection_laplacian(F, G, l, theta, v2fs, L);
}

void get_vertex_connection_laplacian(
                                     const Eigen::MatrixXi & F,
                                     const Eigen::MatrixXi & G,
                                     const Eigen::MatrixXd & l,
                                     const Eigen::MatrixXd & theta,
                                     const Eigen::MatrixXi & v2fs,
                                     Eigen::SparseMatrix<std::complex<double>> & L)
{
    using namespace Eigen;
    using namespace std;
    // Cribbed from igl/cotmatrix_intrinsic

    const int nV = F.maxCoeff()+1;
    L.resize(nV,nV);

    // This is important! it could decrease the comptuation time by a factor of 2
    // Laplacian for a closed 2d manifold mesh will have on average 7 entries per
    // row
    L.reserve(10 * nV);

    // Gather halfedge cotangent weights
    MatrixXd C;
    igl::cotmatrix_entries(l, C);

    // Compute connection. Start by finding halfedge vectors in each vertex
    // MatrixXd theta;
    // MatrixXi v2fs;
    // build_angular_coordinates(F, G, l, theta, v2fs);

    MatrixXcd transportVectorsAlongFaceSide(F.rows(), 3);
    for (int iF = 0; iF < F.rows(); iF++) {
        for (int iS = 0; iS < 3; iS++) {
            Vector2i fs_twin = twin(G, Vector2i{iF, iS});
            double thetaA = theta(iF, iS);
            double thetaB = theta(fs_twin(0), fs_twin(1));
            transportVectorsAlongFaceSide(iF, iS) = -polar(1., thetaB - thetaA);
        }
    }

    vector<Triplet<complex<double>> > T;
    T.reserve(F.rows() * 3 * 4 + nV);

    // Loop over triangles
    for (int iF = 0; iF < F.rows(); iF++) {
        for(int iS = 0; iS < 3; iS++) {
            int src = F(iF, iS);
            int dst = F(iF, (iS+1)%3);
            double w = C(iF, iS);
            complex<double> rot = transportVectorsAlongFaceSide(iF, iS);
            T.push_back(Triplet<complex<double>>(src, src,  w));
            T.push_back(Triplet<complex<double>>(dst, dst,  w));
            T.push_back(Triplet<complex<double>>(src, dst, -w * conj(rot)));
            T.push_back(Triplet<complex<double>>(dst, src, -w * rot));
        }
    }

    double shift = 0; // shift to remove kernel
    for (int iV = 0; iV < nV; iV++) T.push_back(Triplet<complex<double>>(iV, iV, shift));

    L.setFromTriplets(T.begin(), T.end());
}
