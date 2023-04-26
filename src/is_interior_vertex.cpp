#include "is_interior_vertex.h"

bool is_interior_vertex(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v)
{
    bool is_boundary_vertex;
    Eigen::Vector2i fs_v;
    fs_v << v2fs(v,0), v2fs(v,1);
    std::vector<Eigen::Vector2i> one_ring_fs;
    vertex_one_ring_face_sides(G, fs_v, one_ring_fs, is_boundary_vertex);
    if (is_boundary_vertex)
        return false;
    else
        return true;
}