#include "is_v2fs_valid.h"
bool is_v2fs_valid(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & v2fs)
{
    using namespace std;

    int vi, vj;
    Eigen::Vector2i fs;
    int nV = v2fs.rows();

    for (int v=0; v<nV; v++)
    {
        fs << v2fs(v,0) , v2fs(v,1);
        get_face_side_vertices(F,fs,vi,vj);
        if (v != vi) // v2fs must start from v
        {
            cout << "v2fs invalid at v=" << v << endl;
            cout << "v2fs map v to fs=" << fs.transpose() << ", which connects vi=" << vi << "->vj=" << vj << endl;
            return false;
        }
    }
    return true;
}