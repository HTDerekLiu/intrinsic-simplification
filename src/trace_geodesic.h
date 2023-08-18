#ifndef TRACE_GEODESIC
#define TRACE_GEODESIC

#include <Eigen/Core>
#include <vector>

/*
This function traces a geodesic from a point on a mesh intrinsically.

Inputs:
    f_start : integer index of face containing point to start tracing from
    b_start : vector3 of barycentric coordinates of starting point in face f_start
    v_start : vector3 encoding direction to trace, represented in barycentric coordinates
    G       : |F|x6 glue map
    l       : |F|x3 edge lengths for each face side
Outputs:
    f_end   : integer index of face containing ending point
    b_end   : vector3 of barycentric coordinates of endpoing point in face f_end
    path    : if not null, stores the list of edge points hit along the geodesic

    returns true if the tracing terminated successfully, and false if tracing stopped at
    an edge (e.g. due to hitting the boundary)
*/
bool trace_geodesic(
                    int f_start,
                    const Eigen::Vector3d & b_start,
                    const Eigen::Vector3d & v_start,
                    const Eigen::MatrixXi & G,
                    const Eigen::MatrixXd & l,
                    int & f_end,
                    Eigen::Vector3d & b_end,
                    std::vector<std::pair<Eigen::Vector2i, double>>* path = nullptr
                  );
#endif
