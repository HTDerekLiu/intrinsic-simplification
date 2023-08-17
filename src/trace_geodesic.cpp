#include "trace_geodesic.h"

#include <Eigen/Dense>

#include <array>
#include <iostream>

#include <flatten_diamond_mesh.h>

// geodesic tracing helpers
namespace {
    using namespace Eigen;

    const double TRACE_EPS_TIGHT = 1e-12;
    const double TRACE_EPS_LOOSE = 1e-9;
    const bool TRACE_PRINT = false;

    double clamp(double x, double low, double high) {
        return fmin(fmax(x, low), high);
    }

    Vector3d clamp(Vector3d x, double low, double high) {
      return x.cwiseMax(low).cwiseMin(high);
    }

    // given barycentric coordinates which treat side iS as the first side,
    // permute to the canonical ordering
    Vector3d permute_barycentric_to_canonical(Vector3d b, int iS) {
        Vector3d result;
        for (int i = 0; i < 3; i++) result((i+iS)%3) = b(i);
        return result;
    }

    // given barycentric coordinates in the canonical ordering,
    // permute to treat side iS as the first side,
    Vector3d permute_barycentric_from_canonical(Vector3d b, int iS) {
        Vector3d result;
        for (int i = 0; i < 3; i++) result(i) = b((i + iS)%3);
        return result;
    }
}

// (adapted from geometry central's trace_geodesic routines)

// == General routine for tracing a geodesic ray within a face
// If the ray ends within this face, returns true and sets b_end to the
// barycentric coordinates of the endpoint. If not, returns false and sets:
//    s_end   : the index of the side which the ray hits first
//    t_ray   : the time along the ray when the ray hits side s_end
//    t_cross : the time along side s_end when it intersects the ray
//    b_end   : the location of this intersection within face f_start in barycentric coordinates
bool trace_in_face_barycentric(
                               int iF,
                               const Eigen::Vector3d& b_start,
                               const Eigen::Vector3d& v,
                               const std::array<bool, 3> & edge_is_hittable,
                               int & s_end,
                               double & t_ray,
                               double & t_cross,
                               Eigen::Vector3d& b_end)
{
    using namespace Eigen;
    using namespace std;

    if (::TRACE_PRINT) {
        cout << "  general trace in face: " << endl;
        cout << "  face: " << iF << " b_start " << b_start.transpose()
            << " v = " << v.transpose() << endl;
    }

    // Test if the vector ends in the triangle
    b_end = b_start + v;
    if (b_end(0) >= 0 && b_end(1) >= 0 && b_end(2) >= 0) {
        return true;
    }

    // The vector did not end in this triangle. Pick an appropriate point along some edge
    t_ray = numeric_limits<double>::infinity();
    s_end = -1;
    for (int iS = 0; iS < 3; iS++) {
        int iV = (iS + 2) % 3; // vertex opposite side iS

        double t = -b_start(iV) / v(iV); // locate crossing time

        if (!edge_is_hittable[iS] || v(iV) >= 0) continue;

        if (t < t_ray) { // check if we've found the new closest intersection
            t_ray = t;
            s_end = iS;
        }
    }

    t_ray = clamp(t_ray, 0., 1. - ::TRACE_EPS_LOOSE); // clamp to a sane range

    if (s_end < 0) {
        throw runtime_error("trace_in_face_barycentric: no edge intersection found");
    }

    b_end = b_start + t_ray * v; // set b_end to ray-edge intersection
    int v_end = (s_end + 2) % 3; // vertex opposite side s_end
    t_cross = b_end((v_end+2)%3) / (b_end((v_end+1)%3) + b_end((v_end+2)%3));
    return false;
}

Eigen::Matrix3d edge_barycentric_transition_matrix(
                                                   int f_
                                                   );

void trace_geodesic(
                    int f_start,
                    const Eigen::Vector3d & b_start_,
                    const Eigen::Vector3d & v_start_,
                    const Eigen::MatrixXi & G,
                    const Eigen::MatrixXd & l,
                    int & f_end,
                    Eigen::Vector3d& b_end,
                    std::vector<std::pair<Vector2i, double>>* path)
{
    using namespace Eigen;
    using namespace std;
    using namespace global_variables;

    if (::TRACE_PRINT) {
          cout << "\n>>> Trace query (barycentric) from " << f_start << " "
              << b_start_.transpose() << " vec = " << v_start_.transpose() << endl;
    }

    // Early-out if zero
    if (v_start_.squaredNorm() < ::TRACE_EPS_TIGHT) {
        f_end = f_start;
        b_end = b_start_;
    }

    // Make sure the input is sane, i.e. that b_start's components are in the unit interval
    // and sum to 1, and v_start's components sum to zero
    Vector3d b_start = clamp(b_start_, 0, 1);
    b_start /= b_start.sum();
    Vector3d v_start = v_start_ - Vector3d::Constant(v_start_.sum() / 3);

    int f_curr = f_start;
    Vector3d b_curr = b_start;
    Vector3d v_curr = v_start;
    array<bool, 3> hittable{true, true, true};
    int s_curr;
    double t_ray, t_cross;
    while (!trace_in_face_barycentric(f_curr, b_curr, v_curr, hittable,
                                      s_curr, t_ray, t_cross, b_end)) {
        if (path) path->push_back(std::make_pair(Vector2i{f_curr, s_curr}, t_cross));

        //== identify next face
        Vector2i fs_next {G(f_curr, 2 * s_curr), G(f_curr, 2 * s_curr + 1)};
        if (fs_next == GHOST_FACE_SIDE) { // TODO: support barrier edges
          f_end = f_curr;
          return;
        }

        //== subtract off traced component of tangent vector
        v_curr = (1 - t_ray) * v_curr;

        //== update tangent vector
        // first, build map between barycentric coordinate spaces of neighboring triangles
        MatrixXd diamond_layout; MatrixXi ignore;
        // {f, s} goes from diamond_layout.row(0) to diamond_layout.row(1)
        // face f contains diamond_layout.row(2) and its twin contains the final row
        flatten_diamond_mesh(G, l, {f_curr, s_curr}, diamond_layout, ignore);

        Matrix3d V_curr, V_twin; // build position matrices in homogeneous coordinates
        V_curr.col(0) = diamond_layout.row(0);
        V_curr.col(1) = diamond_layout.row(1);
        V_curr.col(2) = diamond_layout.row(2);
        V_curr.row(2) = Vector3d::Ones();
        V_twin.col(0) = diamond_layout.row(1);
        V_twin.col(1) = diamond_layout.row(0);
        V_twin.col(2) = diamond_layout.row(3);
        V_twin.row(2) = Vector3d::Ones();


        if (::TRACE_PRINT) {
            cout << "  converting vector: " << endl;
            cout << "  v (initial) = " << v_curr.transpose() << endl;
        }
        // V_twin^{-1}.V_twin maps from face f to its twin, ordering barycentric coordinates
        // so that side s comes first
        Vector3d position_space = V_curr * permute_barycentric_from_canonical(v_curr, s_curr);
        v_curr = V_twin.colPivHouseholderQr().solve(position_space);
        // project to ensure the vector is in the right direction
        v_curr(2) = fmax(v_curr(2), ::TRACE_EPS_TIGHT);
        // Manual displacement projection to sum to 0 to keep it pointing in the right direction
        double diff = -v_curr.sum();
        if (diff > 0) {
          v_curr(2) += diff;
        } else {
          v_curr += Vector3d{diff, diff, diff} / 3.;
        }
        // reorder barycentric coordinates
        v_curr = permute_barycentric_to_canonical(v_curr, fs_next(1));
        if (::TRACE_PRINT) {
            cout << "  v (transformed) = " << v_curr.transpose() << endl;
        }

        //== update face and barycentric coordinates
        f_curr = fs_next(0);
        int s_curr = fs_next(1);
        b_curr = Vector3d::Zero();
        b_curr((s_curr+0)%3) = t_cross;
        b_curr((s_curr+1)%3) = 1 - t_cross;

        //== set hittable edges
        hittable = array<bool, 3>{true, true, true};
        hittable[s_curr] = false;
    }
    f_end = f_curr;
    b_end = b_curr;
}
