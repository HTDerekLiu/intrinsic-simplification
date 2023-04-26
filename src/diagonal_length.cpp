#include "diagonal_length.h"
double diagonal_length(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs)
{
    using namespace Eigen;
    using namespace std;

    // compute diamond face sides (see .h for notations)
    Vector2i p, q, u, v;
    p = next(next(fs));
    q = next(twin(G, fs));
    u = next(fs);
    v = next(next(twin(G, fs)));

    // compute angles
    double a = opposite_corner_angle(l, p);
    double b = opposite_corner_angle(l, q);

    // compute new edge length
    double lu, lv;
    {
        lu = l(u(0), u(1));
        lv = l(v(0), v(1));
    }
    double diag_edge_length = sqrt(lu*lu + lv*lv - 2*lu*lv*cos(a+b));

    // check positive edge length
    if (diag_edge_length < 0)
    {
        cout << "diagonal_length.cpp got a negative edge length" << endl;
        assert(false);
    }

    return diag_edge_length;
}