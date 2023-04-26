#include "flip_to_delaunay.h"
void flip_to_delaunay(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V)
{
	using namespace std;
    using namespace Eigen;
    
    // A queue of face-sides to test for the Delaunay criterion
    queue<pair<int, int>> to_process;

    // add all face sides to the queue
    int nF = F.rows();
    for (int f=0; f<nF; f++)
        for (int s=0; s<3; s++)
            to_process.push({f, s});

    // An array to keep track whether a fs is in queue
    Matrix<bool,Dynamic,3> is_in_queue;
    is_in_queue.resize(F.rows(), 3);
    is_in_queue.setConstant(true);

    // delaunay flipping
    int num_flips = 0;
    while (to_process.size() > 0)
    {
        // get fs
        pair<int, int> fs_pair = to_process.front(); // access first element
        int f = fs_pair.first;
        int s = fs_pair.second;
        Vector2i fs; fs << f, s;
        to_process.pop(); // remove first eleemnt
        is_in_queue(f, s) = false; // no longer in queue

        // check delaynay
        if (!is_delaunay(G,l,fs) && !is_boundary_face_side(G,fs))
            if (is_diamond_convex(G,l,fs))
            {
                flip_edge(fs,F,G,l,A,v2fs,BC,F2V);
                num_flips += 1;

                // get neighbor fs
                MatrixXi neighbors(2,4);
                neighbors.col(0) << next(fs);
                neighbors.col(1) << next(next(fs));
                neighbors.col(2) << next(twin(G,fs));
                neighbors.col(3) << next(next(twin(G,fs)));

                // enqueue
                for (int ii=0; ii<neighbors.cols(); ii++)
                {
                    int f_n = neighbors(0, ii);
                    int s_n = neighbors(1, ii);
                    if (is_in_queue(f_n, s_n) == false)
                    {
                        to_process.push({ f_n, s_n });
                        is_in_queue(f_n, s_n) = true;
                    }
                }
            }
    }
    // cout << "perform " << num_flips << " edge flips to Delaunay" << endl;
}

void flip_to_delaunay(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs)
{
    Eigen::MatrixXd BC;
    std::vector<std::vector<int>> F2V(F.rows());
    flip_to_delaunay(F,G,l,A,v2fs,BC,F2V);
}