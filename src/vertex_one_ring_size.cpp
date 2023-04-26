#include "vertex_one_ring_size.h"

int vertex_one_ring_size(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G, v2fs, v, fs_list);

    int nF = fs_list.size();
    int one_ring_size = 0;
    
    int vi, vj;
    for (int ii=0; ii<nF; ii++)
    {
        Vector2i fs = fs_list[ii];
        one_ring_size += 1;

        // if encounter boundary face side
        if (is_same_face_side(ccw(G,fs), GHOST_FACE_SIDE))
        {
            get_face_side_vertices(F, next(next(fs)), vi, vj);
            one_ring_size += 1;
        }
    }
    return one_ring_size;
}