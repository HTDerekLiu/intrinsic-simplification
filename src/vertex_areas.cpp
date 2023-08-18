#include "vertex_areas.h"

#include "global_variables.h"

void vertex_areas(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::VectorXd & VA)
{
    using namespace Eigen;

    int nV = F.maxCoeff() + 1;
    int nF = l.rows();

    VectorXd FA;
    face_areas(l, FA);

    VA.resize(nV);
    VA.setZero();
    for (int f=0; f<nF; f++)
    {
        if (F(f, 0) == global_variables::GHOST_INDEX) continue; // skip dead faces
        VA(F(f,0)) += FA(f) / 3;
        VA(F(f,1)) += FA(f) / 3;
        VA(F(f,2)) += FA(f) / 3;
    }
}
