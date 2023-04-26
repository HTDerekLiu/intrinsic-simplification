#include "vertex_one_ring_face_sides.h"

void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs_start,
    std::vector<Eigen::Vector2i> & one_ring_fs)
{
	bool is_boundary_vertex;
	vertex_one_ring_face_sides(G, fs_start, one_ring_fs, is_boundary_vertex);
}

void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs_start,
    std::vector<Eigen::Vector2i> & one_ring_fs,
	bool & is_boundary_vertex)
{
	using namespace std;
    using namespace Eigen;
	using namespace global_variables;

    // track the one-ring face sides
    one_ring_fs.clear();
	one_ring_fs.push_back(fs_start);
	is_boundary_vertex = false; // initialize it to false

	Vector2i fs; // current fs and fs_ccw
	fs << fs_start; // initialize current fs

	while (true)
	{
		fs = ccw(G, fs);
			
		if (is_same_face_side(fs, fs_start))
		{
			// this is an interior vertex and fs_ccw reaches the starting fs
			break;
		}
		if (is_same_face_side(fs, GHOST_FACE_SIDE))
		{
			// fs_ccw hits the boundary, then go clock wise
			is_boundary_vertex = true;
			fs << fs_start;
			Vector2i fs_cw, fs_twin;
			while (true)
			{
				fs_twin = twin(G,fs);
				if (is_same_face_side(fs_twin, GHOST_FACE_SIDE))
				{
					// reaches the right most boundary face side
					break;
				}
				fs = cw(G, fs);
				one_ring_fs.insert(one_ring_fs.begin(), fs);
			}
			break;
		}
		one_ring_fs.push_back(fs);
	}
}

void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
	const int & v,
    std::vector<Eigen::Vector2i> & one_ring_fs)
{
	using namespace global_variables;
	using namespace std;

	bool is_boundary_vertex;
	Eigen::Vector2i fs_start;
	fs_start << v2fs(v,0), v2fs(v,1);
	
	if (!is_same_face_side(fs_start, GHOST_FACE_SIDE)) // if fs_start is valid
		vertex_one_ring_face_sides(G, fs_start, one_ring_fs, is_boundary_vertex);
	else // fs_start is invalid
		one_ring_fs.clear(); // make sure one_ring_fs is empty
}

void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
	const int & v,
    std::vector<Eigen::Vector2i> & one_ring_fs,
	bool & is_boundary_vertex)
{
	using namespace global_variables;
	using namespace std;

	Eigen::Vector2i fs_start;
	fs_start << v2fs(v,0), v2fs(v,1);
	
	if (!is_same_face_side(fs_start, GHOST_FACE_SIDE)) // if fs_start is valid
		vertex_one_ring_face_sides(G, fs_start, one_ring_fs, is_boundary_vertex);
	else // fs_start is invalid
		one_ring_fs.clear(); // make sure one_ring_fs is empty
}