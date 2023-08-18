#ifndef INSERT_EAR_VERTEX
#define INSERT_EAR_VERTEX

#include <Eigen/Core>

/*
    This function inserts an ear vertex into a boundary face side

Inputs:
    fs    face side in which to insert the new vertex. Should be the interior
            face side of a boundary edge
    t:    barycentric coordinate along fs to place the new vertex at
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
int insert_ear_vertex(
                      const Eigen::Vector2i & fs,
                      double t,
                      Eigen::MatrixXi & F,
                      Eigen::MatrixXi & G,
                      Eigen::MatrixXd & l,
                      Eigen::MatrixXd & A,
                      Eigen::MatrixXi & v2fs,
                      Eigen::MatrixXd & BC,
                      std::vector<std::vector<int>> & F2V);

#endif
