#include "is_boundary_face_side.h"

bool is_boundary_face_side(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs)
{
    using namespace Eigen;
	using namespace global_variables;

	Vector2i fs_twin = twin(G, fs);
	return is_same_face_side(fs_twin, GHOST_FACE_SIDE);
}