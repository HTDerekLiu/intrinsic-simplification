#include "is_ear_vertex.h"

bool is_ear_vertex(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v)
{
    using namespace Eigen;
    using namespace std;

    if (is_interior_vertex(G,v2fs,v))
        return false;
    else // is boundary vertex
    {
        vector<Vector2i> fs_list;
        vertex_one_ring_face_sides(G, v2fs, v, fs_list);
        if (fs_list.size() == 1)
            return true;
        else
            return false;
    }
}