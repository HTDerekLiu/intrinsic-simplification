#include "get_angular_coordinate_range.h"
double get_angular_coordinate_range(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & v2fs,
    const int & v)
{
    using namespace std;
    using namespace Eigen;

    if (is_boundary_vertex(G,v2fs,v))
    {
        vector<Vector2i> fs_list;
        vertex_one_ring_face_sides(G, v2fs, v, fs_list);
        if (fs_list.size() > 1) // regular boundary vertex
        {
            double angle_sum = 0.0;
            for (int ii=0; ii<fs_list.size(); ii++)
                angle_sum += opposite_corner_angle(l, next(fs_list[ii]));
            return angle_sum / (angle_sum+M_PI) * 2 * M_PI;
            
        }
        else // ear vertex
            return M_PI;
    }
    else // interior vertex
        return 2 * M_PI;
    
}