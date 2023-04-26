#include "face_sides_to_curve_network.h"

void face_sides_to_curve_network(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const std::vector<Eigen::Vector2i> & fs_list,
    Eigen::MatrixXd & nodes,
    Eigen::MatrixXi & edges)
{
    using namespace Eigen;

    int num_fs = fs_list.size();
 	nodes.resize(num_fs*2, 3);
	edges.resize(num_fs, 2);

	for (int ii=0; ii<num_fs; ii++)
	{
		Vector2i fs = fs_list[ii];
		int vi, vj;
		get_face_side_vertices(F,fs,vi,vj);

		nodes.row(2*ii) = V.row(vi);
		nodes.row(2*ii+1) = V.row(vj);
		edges.row(ii) << 2*ii, 2*ii+1;
	}

}