#include "angle_sum.h"

void angle_sum(
    const Eigen::MatrixXi &F, 
    const Eigen::MatrixXd &l, 
    Eigen::VectorXd &ang_sum)
{
    int nV = F.maxCoeff() + 1;
    angle_sum(nV,F,l,ang_sum);
}

void angle_sum(
    const int & nV,
    const Eigen::MatrixXi &F, 
    const Eigen::MatrixXd &l, 
    Eigen::VectorXd &ang_sum)
{
    using namespace Eigen;
    using namespace global_variables;

    ang_sum.resize(nV);
    // ang_sum.setZero();
    ang_sum.setConstant(DOUBLE_INF);

    for (int f=0; f<F.rows(); f++){
        for (int s=0; s<3; s++){
            int v = F(f,s);
            if (v != GHOST_INDEX){
                Vector2i fs = {f, s};
                Vector2i fs_opp = next(next(fs));
                double theta = opposite_corner_angle(l, fs_opp);
                if (ang_sum(v) == DOUBLE_INF){
                    ang_sum(v) = 0.0;
                }
                ang_sum(v) += theta;
            }
        }
    }
}