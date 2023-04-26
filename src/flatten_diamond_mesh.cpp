#include "flatten_diamond_mesh.h"

void flatten_diamond_mesh(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs,
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & F)
{
    using namespace Eigen;

    // flatten the diamond
    U.resize(4,2);
    U.setZero();

    double scale = 1.0 / l(fs(0), fs(1)); // scale edges to avoid numerical issue
    U.row(0) << 0.0, 0.0;
    U.row(1) << 0.0, 1.0;

    // collect useful face sides
    Vector2i fs_next = next(fs);
    Vector2i fs_next_next = next(fs_next);
    Vector2i fs_twin = twin(G, fs);
    Vector2i fs_twin_next = next(fs_twin);
    Vector2i fs_twin_next_next = next(fs_twin_next);

    // compute U.row(2)
    {
        double ang2 = opposite_corner_angle(l, fs_next);
        ang2 = M_PI/2.0 - ang2;

        double e02 = l(fs_next_next(0), fs_next_next(1));
        U.row(2) << -cos(ang2)*e02*scale, sin(ang2)*e02*scale;
    }

    // compute U.row(3)
    {
        double ang3 = opposite_corner_angle(l, fs_twin_next_next);
        ang3 = M_PI/2.0 - ang3;

        double e03 = l(fs_twin_next(0), fs_twin_next(1));
        U.row(3) << cos(ang3)*e03*scale, sin(ang3)*e03*scale;
    }

    // compute F
    F.resize(2,3);
    {
        Vector3i tmp, tmp_rolled;

        tmp << 0, 1, 2;
        roll1d(tmp, fs(1), tmp_rolled);
        F.row(0) << tmp_rolled.transpose();

        tmp << 1, 0, 3;
        roll1d(tmp, fs_twin(1), tmp_rolled);
        F.row(1) << tmp_rolled.transpose();
    }
}

void flatten_diamond_mesh(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs,
    const bool & is_ccw, // is "ccw" or "cw" flip direction
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & F_flip)
{
    using namespace Eigen;
    flatten_diamond_mesh(G, l, fs, U, F);

    // create face list for the flipped triangulation
    F_flip.resize(2,3);
    Vector2i fs_twin = twin(G, fs);
    if (is_ccw) // flip counter clock wise
    {
        Vector3i tmp, tmp_rolled;

        tmp << 3, 2, 0;
        roll1d(tmp, fs(1), tmp_rolled);
        F_flip.row(0) << tmp_rolled.transpose();

        tmp << 2, 3, 1;
        roll1d(tmp, fs_twin(1), tmp_rolled);
        F_flip.row(1) << tmp_rolled.transpose();
    }
    else // flip clock wise
    {
        Vector3i tmp, tmp_rolled;

        tmp << 2, 3, 1;
        roll1d(tmp, fs(1), tmp_rolled);
        F_flip.row(0) << tmp_rolled.transpose();

        tmp << 3, 2, 0;
        roll1d(tmp, fs_twin(1), tmp_rolled);
        F_flip.row(1) << tmp_rolled.transpose();
    }
}