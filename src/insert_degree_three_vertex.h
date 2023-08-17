#ifndef INSERT_DEGREE_THREE_VERTEX
#define INSERT_DEGREE_THREE_VERTEX

#include <Eigen/Core>

/*
    This function inserts a degree three vertex into a face

Inputs:
    f:    face in which to insert the new vertex
    b:    barycentric coordinates for new vertex in face iF
    F:    |F|x3 vertex-face adjacency list
    G:    |F|x6 glue map
    l:    |F|x3 edge lengths for each face side
    A:    |F|x3 angular coordinate for each face side
    v2fs: |V|x2 where v2fs.row(i) returns a face side for vertex i
    BC:   |BC|x3 array of barycentric coordinates whose corresponding faces are stored
            in F2V implicitly
    F2V:  |F| length list of lists, where F2V[f] gives you a list of indices in BC.
            For example, if F2V[f] = [v], then BC[v,:] corresponds to the barycentric
            coordinates in F[f,:]

Outputs:
    F, G, l, A, v2fs, BC, F2V, and T are changed in place
    returns the index of the inserted vertex
*/
int insert_degree_three_vertex(
    int f,
    const Eigen::Vector3d & b,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V);

#endif
