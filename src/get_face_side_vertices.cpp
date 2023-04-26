#include "get_face_side_vertices.h"

void get_face_side_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::Vector2i & fs,
    int & vi,
    int & vj)
{
    int f = fs(0);
    int s = fs(1);
    vi = F(f, s);
    vj = F(f, (s+1)%3);
};