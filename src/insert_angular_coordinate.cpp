#include "insert_angular_coordinate.h"

void insert_angular_coordinate(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs,
    Eigen::MatrixXd & A)
{
    using namespace Eigen;
    using namespace std;
    using namespace global_variables;

    // insert the angular coordinate for s0 and twin(s0)
	Vector2i fs_next = next(fs);
	Vector2i fs_twin_next_next = next(next(twin(G, fs)));

	double left_angle = opposite_corner_angle(l, fs_next);
	double right_angle = opposite_corner_angle(l, fs_twin_next_next);

	double ratio = right_angle / (right_angle + left_angle);

	Vector2i fs_cw; 
	double angle_cw;
	{
		fs_cw = cw(G, fs);
		angle_cw = get_face_side_angular_coordinate(A, fs_cw);
	}

	Vector2i fs_ccw;
	double angle_ccw;
	{
		fs_ccw = ccw(G, fs);
		
		// extract angle_ccw, but fs_ccw could be outside the mesh
		if (is_same_face_side(fs_ccw, GHOST_FACE_SIDE))
		{
			// if fs_ccw is outside the mesh
			vector<Vector2i> fs_list;
			vertex_one_ring_face_sides(G, fs, fs_list);
			
			VectorXd angles(fs_list.size());
			for (auto it = std::begin(fs_list); it != std::end(fs_list); ++it)
			{
				Vector2i fs_next = next(*it);

				double angle = opposite_corner_angle(l, fs_next);

				int idx = it - fs_list.begin();
				angles[idx] = angle;
			}

			angle_ccw = angles.sum() / (angles.sum()+M_PI) * M_PI*2;
		}
		else
		{
			// if fs_ccw is inside the mesh
			angle_ccw = A(fs_ccw(0), fs_ccw(1));
		}
		
		if (angle_ccw < angle_cw)
        	angle_ccw += 2*M_PI;
		
		double angle_diff = angle_ccw - angle_cw;
		double angular_coord_fs = remainder((angle_cw + ratio * angle_diff), (2*M_PI));
		
		// assign fs to angular coordinate
		if (!is_same_face_side(fs, GHOST_FACE_SIDE))
			A(fs(0), fs(1)) = angular_coord_fs;
	}
}