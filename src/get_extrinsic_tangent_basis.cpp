#include "get_extrinsic_tangent_basis.h"

void get_extrinsic_tangent_basis(
                                 const Eigen::MatrixXd & V,
                                 const Eigen::MatrixXi & F,
                                 const Eigen::MatrixXd & A,
                                 const Eigen::MatrixXi & v2fs,
                                 Eigen::MatrixXd & basisX,
                                 Eigen::MatrixXd & basisY)
{

    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, N);

    basisX.resize(V.rows(), 3);
    basisY.resize(V.rows(), 3);

    for (int iV = 0; iV < V.rows(); iV++) {
        int iF = v2fs(iV, 0);
        int iS = v2fs(iV, 1);
        Eigen::Vector3d vN = N.row(iV);
        Eigen::Vector3d vEdge = V.row(F(iF, (iS + 1)%3)) - V.row(F(iF, iS));

        // first compute a basis where v.halfedge() is the x-axis
        Eigen::Vector3d vX = vEdge - vN.dot(vEdge) * vEdge;
        vX.normalize();
        Eigen::Vector3d vY = vN.cross(vX);

        // now take angular coordinate into account.
        // really v.halfedge() should be at an angle of A(iF, iS) from the x axis
        double theta = A(iF, iS);
        basisX.row(iV) = cos(theta) * vX - sin(theta) * vY;
        basisY.row(iV) = sin(theta) * vX + cos(theta) * vY;
    }
}
