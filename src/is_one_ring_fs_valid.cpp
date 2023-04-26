#include "is_one_ring_fs_valid.h"

bool is_one_ring_fs_valid(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v)
{
    using namespace std;

    std::vector<Eigen::Vector2i> one_ring_fs;
    vertex_one_ring_face_sides(G, v2fs, v, one_ring_fs);

    bool is_valid = true;
    for (int ii=0; ii<one_ring_fs.size(); ii++)
    {
        int vi, vj;
        get_face_side_vertices(F,one_ring_fs[ii],vi, vj);
        if (vi != v)
        {
            is_valid = false;
            break;
        }
    }

    if (!is_valid)
    {
        cout << "=======\n";
        cout << "one-ring face sides of vertex " << v << " is invalid\n";
        cout << "all one-ring face sides: \n";
        for (int ii=0; ii<one_ring_fs.size(); ii++)
        {
            int vi, vj;
            get_face_side_vertices(F,one_ring_fs[ii],vi, vj);

            int f = one_ring_fs[ii](0);
            int s = one_ring_fs[ii](1);
            cout << "(" << f << "," << s << "): " << vi << " -> " << vj << endl;
        }
        cout << "all one-ring faces: \n";
        for (int ii=0; ii<one_ring_fs.size(); ii++)
        {
            int f = one_ring_fs[ii](0);
            cout << "F " << f << ": " << F(f,0) << "," << F(f,1) << "," << F(f,2) << endl;
        }
        cout << "=======\n";
    }
    return is_valid;
}