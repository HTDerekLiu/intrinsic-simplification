#include "is_diamond_convex.h"

bool is_diamond_convex(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs)
{
    using namespace Eigen;
	using namespace std;
    using namespace global_variables;

	Vector2i fs_next = next(fs);
	Vector2i fs_next_next = next(fs_next);
	Vector2i fs_twin = twin(G, fs);
	Vector2i fs_twin_next = next(fs_twin);
	Vector2i fs_twin_next_next = next(fs_twin_next);

	double angle_i_jk = opposite_corner_angle(l, fs_next);
	double angle_i_lj = opposite_corner_angle(l, fs_twin_next_next);

	double angle_j_ki = opposite_corner_angle(l, fs_next_next);
	double angle_j_il = opposite_corner_angle(l, fs_twin_next);

	double angle_i = angle_i_jk + angle_i_lj;
	double angle_j = angle_j_ki + angle_j_il;
	return (angle_i<(M_PI+EPS) && angle_j<(M_PI+EPS));
}