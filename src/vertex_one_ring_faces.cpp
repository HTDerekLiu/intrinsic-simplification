#include "vertex_one_ring_faces.h"

void vertex_one_ring_faces(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_faces)
{
    using namespace std;
    using namespace Eigen; 

    // create local face sides
    vector<Vector2i> fs_list;
    {
        vertex_one_ring_unique_face_sides(G, Vector2i(v2fs(v, 0), v2fs(v, 1)), fs_list);
    }

    // create local faces
    int nF = fs_list.size();
    one_ring_faces.resize(nF);
    {
        for (int ii=0; ii<nF; ii++)
            one_ring_faces(ii) = fs_list[ii](0);
    }
}