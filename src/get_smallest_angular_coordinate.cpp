#include "get_smallest_angular_coordinate.h"

Eigen::Vector2i get_smallest_angular_coordinate(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & A,
    const Eigen::Vector2i & fs)
{
    using namespace std;
    using namespace Eigen;
    
    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G, fs, fs_list);

    double a_min = A(fs(0), fs(1));
    Vector2i fs_min;
    fs_min << fs;

    for (auto it = std::begin(fs_list); it != std::end(fs_list); ++it)
    {
        Vector2i fs_cur = *it; 
        double a = A(fs_cur(0), fs_cur(1));
        if (a < a_min)
        {
            a_min = a;
            fs_min << fs_cur;
        }
    }

    return fs_min;
}