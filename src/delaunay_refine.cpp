#include "delaunay_refine.h"

#include "next.h"
#include "twin.h"
#include "opposite_corner_angle.h"
#include "face_circumradius.h"
#include "is_delaunay.h"
#include "is_boundary_face_side.h"
#include "is_diamond_convex.h"
#include "is_interior_vertex.h"
#include "always_flip_to_degree_three.h"
#include "always_flip_to_ear.h"
#include "flip_edge.h"
#include "remove_degree_three_vertex.h"
#include "remove_ear_vertex.h"
#include "is_ear_vertex.h"
#include "trace_geodesic.h"
#include "insert_degree_three_vertex.h"
#include "insert_ear_vertex.h"

#include <queue>
#include <array>
#include <vector>
#include <unordered_map>

// find all vertices within a distance ball_radius of vertex iV in the
// edge graph of the given mesh
std::unordered_map<int, double> vertex_dijkstra_ball(Eigen::MatrixXi & F,
                                                        Eigen::MatrixXi & G,
                                                        Eigen::MatrixXd & l,
                                                        Eigen::MatrixXi & v2fs,
                                                        int iV,
                                                        double ball_radius)
{
    using namespace Eigen;
    using namespace std;

    // Search state: distances discovered so far
    unordered_map<int, double> shortest_dist;

    // Search state: visible neighbors eligible to expand to
    using w_vertex = tuple<double, int>; // weighted vertex
    priority_queue<w_vertex, vector<w_vertex>, greater<w_vertex>> to_process;

    // Add initial source
    to_process.emplace(0., iV);

    while (!to_process.empty()) {

        // Get the next closest neighbor off the queue
        double curr_dist = get<0>(to_process.top());
        int    curr_vert = get<1>(to_process.top());
        to_process.pop();

        // skips stale entries
        if (shortest_dist.find(iV) != shortest_dist.end()) continue;

        shortest_dist[curr_vert] = curr_dist; // lock in distance

        // propagate distance to all of its neighbors
        vector<Vector2i> outgoing_face_sides;
        vertex_one_ring_face_sides(G, v2fs, curr_vert, outgoing_face_sides);
        for (Vector2i fs : outgoing_face_sides) {
            double target_dist = curr_dist + l(fs(0), fs(1));
            int    target_vert = F(fs(0), (fs(1) + 1) %3);
            if (target_dist <= ball_radius
                && shortest_dist.find(target_vert) == shortest_dist.end()) {

                to_process.emplace(target_dist, target_vert);
            }
        }
    }

    return shortest_dist;
}

void delaunay_refine(
                     Eigen::MatrixXi & F,
                     Eigen::MatrixXi & G,
                     Eigen::MatrixXd & l,
                     Eigen::MatrixXd & A,
                     Eigen::MatrixXi & v2fs,
                     Eigen::MatrixXd & BC,
                     std::vector<std::vector<int>> & F2V,
                     double angle_threshold_degrees,
                     double circumradius_threshold,
                     int max_iterations)
{
    using namespace Eigen;
    using namespace std;
    using namespace global_variables;
    using VectorXb = Matrix<bool, Dynamic, 1>;
    using MatrixXb = Matrix<bool, Dynamic, Dynamic>;

    int n_orig_vertices = F.maxCoeff() + 1;
    auto is_inserted_vertex = [&](int iV) -> bool {
        return iV >= n_orig_vertices;
    };

    // TODO: do we need this?
    // initialize transport cost
    MatrixXd T = MatrixXd::Zero(n_orig_vertices, 9);

    // borrowed from geometry central
    // Relationship between angles and circumradius-to-edge
    double angle_threshold_radians = angle_threshold_degrees * M_PI / 180.;

    // store whether an edge is flippable
    // TODO: mark whether boundary edges are flippable, and update following
    //       mesh mutations?
    MatrixXb is_fixed(l.rows(), l.cols());
    is_fixed.fill(false);

    // TODO: do we need this?
    // store whether a face still exists
    // (although we mostly refine the mesh, faces can be deleted when we remove
    //  inserted vertices in a ball following an edge split)
    // VectorXb is_dead(F.rows());
    // is_dead.fill(false);

    // shorthand for twin function using glue map G
    auto twin = [&](Vector2i fs) -> Vector2i {
      return ::twin(G, fs);
    };

    // compute area by Heron's rule
    auto face_area = [&](size_t iF) -> double {
        double a = l(iF, 0);
        double b = l(iF, 1);
        double c = l(iF, 2);
        double s = (a + b + c) / 2.0;
        double squared_area = s * (s - a) * (s - b) * (s - c);
        squared_area = fmax(0., squared_area); // clamp squared area to be positive
        return sqrt(squared_area);
    };

    // test if a face violates the circumradius ratio condition
    auto needs_refinement = [&](int iF) {
        size_t n_needle = 0;
        for (size_t iS = 0; iS < 3; iS++) {
            // compute vertex angle sum by looping over outgoing face sides
            double angle_sum = 0;
            Vector2i fs_start{iF, iS};
            Vector2i fs_curr = fs_start;
            // cout << "Loop 1 begin" << endl;
            do {
                angle_sum += opposite_corner_angle(l, next(fs_curr));
                fs_curr = twin(next(next(fs_curr)));
            } while (fs_curr != fs_start && !is_fixed(fs_curr(0), fs_curr(1)));
            // cout << "Loop 1 end" << endl;

            // traverse other direction if we stopped due to hitting a fixed edge
            if (fs_curr != fs_start && !is_fixed(fs_start(0), fs_start(1))) {
                // cout << "Loop 2 begin" << endl;
                fs_curr = fs_start;
                do {
                    fs_curr = next(twin(fs_curr));
                    angle_sum += opposite_corner_angle(l, next(fs_curr));
                } while (!is_fixed(fs_curr(0), fs_curr(1)));
                // cout << "Loop 2 end" << endl;
            }

            if (angle_sum < M_PI / 3.) n_needle++;
        }
        // cout << "face check end" << endl;
        if (n_needle == 1) return false;

        // TODO: can we implement this check?
        // Face inputFace = getParentFace(f);
        // if (inputFace != Face()) {
        //     inputGeom.requireVertexAngleSums();
        //     for (Vertex v : inputFace.adjacentVertices()) {
        //         if (inputGeom.vertexAngleSums[v] < M_PI / 3.) {
        //             inputGeom.unrequireVertexAngleSums();
        //             return false;
        //         }
        //     }
        //     inputGeom.unrequireVertexAngleSums();
        // }

        double fc = face_circumradius(l, iF);
        // double fl = shortest_edge(iF);

        bool needs_refinement_length = fc > circumradius_threshold;

        // Explicit check allows us to skip degree one vertices
        // (can't make those angles smaller!)
        bool needs_refinement_angle = false;
        for (size_t iS = 0; iS < 3; iS++) {
            Vector2i fs{iF, iS};
            double base_angle = opposite_corner_angle(l, next(fs));
            if (base_angle < angle_threshold_radians) {
                // If it's already a degree one vertex, nothing we can do here
                bool is_degree_one_vertex = (twin(fs) == next(next(fs)));
                if (is_degree_one_vertex) {
                    continue;
                }

                // If it's a fixed corner, can't make it smaller
                if (is_fixed(iF, iS) && is_fixed(iF, (iS+2)%3)) {
                    continue;
                }

                needs_refinement_angle = true;
            }
        }

        return needs_refinement_length || needs_refinement_angle;
    };

    // Manages a check at the bottom to avoid infinite-looping when
    // numerical baddness happens
    int recheck_count = 0;
    const int MAX_RECHECK_COUNT = 5;

    // Track statistics
    size_t n_flips = 0;
    size_t n_insertions = 0;

    // Initialize queue of (possibly) non-delaunay edges by pushing all
    // edges onto the queue
    queue<Vector2i> delaunay_check_queue;
    MatrixXb in_delaunay_queue(F.rows(), 3);
    in_delaunay_queue.setConstant(true);

    for (size_t iF = 0; iF < F.rows(); iF++) {
        for (size_t iS = 0; iS < 3; iS++) {
            delaunay_check_queue.push({iF, iS});
        }
    }

    // check if fs or its twin is in delaunay queue
    auto fs_in_delaunay_queue = [&](Vector2i fs) -> bool {
        Vector2i fsT = twin(fs);
        return in_delaunay_queue(fs(0), fs(1))
          || (fsT != GHOST_FACE_SIDE && in_delaunay_queue(fsT(0), fsT(1)));
    };

    // Return a weight to use for sorting PQ. Usually sorts by biggest area,
    // but also puts faces on the boundary first with weight inf
    auto area_weight = [&](int iF) -> double {
      for (size_t iS = 0; iS < 3; iS++) {
          if (is_fixed(iF, iS)) return numeric_limits<double>::infinity();
      }
      return face_area(iF);
    };

    // Define a queue of (possibly) circumradius-violating faces, processing the
    // largest faces first (good heuristic)
    typedef tuple<double, double, int> a_face;
    priority_queue<a_face, vector<a_face>, less<a_face>> circumradius_check_queue;

    auto enqueue_face_circumradius = [&](int iF) {
        circumradius_check_queue.push(make_tuple(area_weight(iF), face_area(iF), iF));
    };

    // Function to check the neighbors of an edge for further processing after a flip.
    // TODO: somehow check neighbors after flipping to degree 3 for vertex removal
    auto check_neighbors_after_flip = [&](Vector2i fs) {
        // cout << "  flipped face side " << fs << endl;
        n_flips++;

        // Add neighboring faces, which might violate circumradius constraint
        array<int, 2> neighboring_faces{fs(0), twin(fs)(0)};
        for (int nF : neighboring_faces) {
            if (needs_refinement(nF)) enqueue_face_circumradius(nF);
        }

        // Add neighbors to delaunay queue, as they may need flipping now
        Vector2i fsN = next(fs);
        Vector2i fsT = twin(fs);
        Vector2i fsTN = next(fsT);
        array<Vector2i, 4> neigh_edges {fsN, next(fsN), fsTN, next(fsTN)};
        for (Vector2i fsN : neigh_edges) {
            if (!fs_in_delaunay_queue(fsN)) {
                delaunay_check_queue.push(fsN);
                in_delaunay_queue(fsN(0), fsN(1)) = true;
            }
        }
    };

    // Flip the triangulation back to being Delaunay, using the
    // carefully-maintained queue of possibly non-Delaunay edges.
    auto flip_to_delaunay_from_queue = [&]() {
        while (!delaunay_check_queue.empty()) {

            // Get the top element from the queue of possibily non-Delaunay face sides
            Vector2i fs = delaunay_check_queue.front();
            delaunay_check_queue.pop();
            // if (is_dead(fs(0))) continue; // TODO: do we need this?
            if (!in_delaunay_queue(fs(0), fs(1)))
                throw runtime_error("popped face side which is not in queue");
            in_delaunay_queue(fs(0), fs(1)) = false;

            // flip if flippable and not delaunay
            if (!is_delaunay(G, l, fs) && !is_boundary_face_side(G, fs)) {
                if (is_diamond_convex(G,l,fs)) {
                    flip_edge(fs, F, G, l, A, v2fs, BC, F2V);
                    check_neighbors_after_flip(fs);
                }
            }
        }
    };

    // flip to Delaunay and then initialize queue
    flip_to_delaunay_from_queue();

    for (size_t iF = 0; iF < F.rows(); iF++) {
        if (needs_refinement(iF)) enqueue_face_circumradius(iF);
    }

    auto remove_inserted_vertex = [&](int v) -> int {
        // TODO: grab resulting face to return
        if (is_interior_vertex(G,v2fs,v)) {
            always_flip_to_degree_three(v,F,G,l,A,v2fs,BC,F2V);
            remove_degree_three_vertex(v,F,G,l,A,v2fs,BC,F2V,T);
        } else if (is_ear_vertex(G,v2fs,v)) {
            remove_ear_vertex(v,F,G,l,A,v2fs,BC,F2V,T);
        } else {// regular boundary vertex
            always_flip_to_ear(v,F,G,l,A,v2fs,BC,F2V);
            remove_ear_vertex(v,F,G,l,A,v2fs,BC,F2V,T);
        }
        return -1; // TODO: return new face
    };

    // Subroutine to be invoked to delete previously-inserted vertices whenever
    // refinment splits an edge
    auto delete_nearby_vertices = [&](Vector2i fs1, Vector2i fs2) {
        // radius of the diametral ball
        double ball_radius = max(l(fs1(0), fs1(1)), l(fs2(0), fs2(1)));
        size_t new_v = F(fs1(0), fs1(1));

        // Flip to Delaunay, to ensure that the Dijkstra search below actually
        // has a stretch factor of 2
        flip_to_delaunay_from_queue();

        // Find all vertices within range.
        // Most properly, this should probably be a polyhedral geodesic ball search,
        // but that creates a dependence on polyhedral shortest paths which is bad
        // for performance and robustness
        //
        // Fortunately, on a Delaunay triangulation, the Dijkstra distance is at
        // most 2x the geodesic distance (see Intrinsic Triangulations Course, the
        // underlying reference is Ge Xia 2013. "The Stretch Factor of the Delaunay
        // Triangulation Is Less than 1.998").
        // So instead, we delete all previously-inserted vertices within 2x the
        // Dijkstra radius instead. This may delete some extra verts, but that does
        // not affect convergence.
        unordered_map<int, double> nearby_verts
          = vertex_dijkstra_ball(F, G, l, v2fs, new_v, 2. * ball_radius);

        // remove inserted vertices
        for (auto p : nearby_verts) {
            int iV = p.first;
            if (is_inserted_vertex(iV)) {
                int iF = remove_inserted_vertex(iV);

                if (iF >= 0) {
                    // Add adjacent edges for Delaunay check
                    for (size_t iS = 0; iS < 3; iS++) {
                        if (!fs_in_delaunay_queue({iF, iS})) {
                            delaunay_check_queue.push({iF, iS});
                            in_delaunay_queue(iF, iS) = true;
                        }
                    }

                    // Add face for circumcircle check
                    if (needs_refinement(iF))
                        enqueue_face_circumradius(iF);
                }
            }
        }
    };

    // Insert the a vertex located at the intrinsic circumcenter of the input face.
    // Mesh must be Delaunay. Returns the index of the inserted vertex (or -1 in
    // case of error)
    auto insert_circumcenter = [&](int iF) -> int {
        double a = l(iF, 1);
        double b = l(iF, 2);
        double c = l(iF, 0);
        double a2 = a * a;
        double b2 = b * b;
        double c2 = c * c;
        Vector3d circumcenter_loc{a2 * (b2 + c2 - a2),
                                  b2 * (c2 + a2 - b2),
                                  c2 * (a2 + b2 - c2)};
        circumcenter_loc /= circumcenter_loc.sum();

        // Trace from the barycenter (have to trace from somewhere)
        Vector3d barycenter = Vector3d::Constant(1. / 3.);
        Vector3d vec_to_circumcenter = circumcenter_loc - barycenter;

        // === Trace the ray to find the location of the new point

        // face and barycentric coordinates for circumcenter
        int f_c; Vector3d b_c;
        bool end_in_face = trace_geodesic(iF, barycenter, vec_to_circumcenter, F, G,
                                          l, f_c, b_c);

        if (end_in_face) {
            // ===  Add the new vertex
            return insert_degree_three_vertex(f_c, b_c, F, G, l, A, v2fs, BC, F2V);
        } else {
            // If the circumcenter is blocked by an edge, insert the midpoint of
            // that edge instead (which is needed for Chew's 2nd algo).

            // identify hit edge
            int iS = (is_fixed(f_c, 0) && b_c(2) < 1e-8) ? 0
                   : (is_fixed(f_c, 1) && b_c(0) < 1e-8) ? 1
                   : (is_fixed(f_c, 2) && b_c(1) < 1e-8) ? 2
                   : -1;

            if (iS < 0) {
                throw runtime_error("geodesic trace terminated at non-fixed edge?");
            }

            return insert_ear_vertex({f_c, iS}, 0.5,
                                      F, G, l, A, v2fs, BC, F2V);
        }
        /* TODO
        if (newPositionOnIntrinsic.type == SurfacePointType::Edge) {
        newPositionOnIntrinsic.tEdge = 0.5;
        }
        */

        // cout << " ... inserting vertex in face " << f_circumcenter << " at position "
        //      << b_circircumcenter.transpose() << endl;

    };

    // === Outer iteration: flip and insert until we have a mesh that satisfies both angle
    //                      and circumradius goals
    // static int iInsert = 0;
    do {

        // == First, flip to delaunay
        // cout << "flipping... queue has size " << delaunay_check_queue.size() << endl;
        flip_to_delaunay_from_queue();
        // cout << "done flipping... queue has size " << delaunay_check_queue.size() << endl;

        // check Delaunay condition
        if (false) {
            for (int iF = 0; iF < F.rows(); iF++) {
                for (int iS = 0; iS < 3; iS++) {
                    if (fs_in_delaunay_queue({iF, iS})) {
                      Vector2i fs{iF, iS};
                      Vector2i fsT = twin(fs);
                      cout << "fs: "<< fs.transpose() << " | twin: " << fsT.transpose() << endl;
                      cout << boolalpha << "fs in queue: " << in_delaunay_queue(iF, iS)
                          << " | twin in queue: " << (fsT != GHOST_FACE_SIDE && in_delaunay_queue(fsT(0), fsT(1))) << endl;
                      throw runtime_error("face side left in queue");
                    }
                    if (!is_delaunay(G, l, {iF, iS})) {
                        cout << "twin: " << twin({iF, iS}).transpose() << endl;
                        throw runtime_error("delaunay flipping failed on face side {" + to_string(iF) + ", " + to_string(iS) + "}");
                    }
                }
            }
        }

        // == Second, insert one circumcenter

        // If we've already inserted the max number of points, call it a day
        if (max_iterations >= 0 && n_insertions == max_iterations) {
            break;
        }

        // Try to insert just one circumcenter
        if (!circumradius_check_queue.empty()) {
            // Get the biggest face
            int iF = get<2>(circumradius_check_queue.top());
            double A = get<1>(circumradius_check_queue.top());
            circumradius_check_queue.pop();
            // if (f.isDead()) continue; // TODO: do we need this?

            // Two things might have changed that would cause us to skip this entry:
            //   - If the area has changed since this face was inserted in to the queue, skip it.
            //     Note that we don't need to re-add it, because it must have been placed in the
            //     queue when its area was changed
            //   - This face might have been flipped to no longer violate constraint
            if (A == face_area(iF) && needs_refinement(iF)) {
                // cout << "inserting circumcenter " << iInsert++ << endl;
                int new_v = insert_circumcenter(iF);
                in_delaunay_queue.conservativeResize(F.rows(), 3); // resize buffer to fit two new faces!
                in_delaunay_queue(F.rows()-2, 0) = false;
                in_delaunay_queue(F.rows()-2, 1) = false;
                in_delaunay_queue(F.rows()-2, 2) = false;
                in_delaunay_queue(F.rows()-1, 0) = false;
                in_delaunay_queue(F.rows()-1, 1) = false;
                in_delaunay_queue(F.rows()-1, 2) = false;
                if (new_v < 0) {
                    // vertex insertion failed (probably due to a tracing error)
                    continue;
                }
                n_insertions++;

                // Mark everything in the 1-ring as possibly non-Delaunay and possibly violating
                // the circumradius constraint
                vector<Vector2i> outgoing_face_sides;
                vertex_one_ring_face_sides(G, v2fs, new_v, outgoing_face_sides);
                for (Vector2i fs : outgoing_face_sides) {
                    int iF = fs(0);
                    // check circumradius constraint
                    if (needs_refinement(iF)) {
                        circumradius_check_queue.push(make_tuple(area_weight(iF), face_area(iF), iF));
                    }

                    // Check delaunay constraint
                    for (int iS = 0; iS < 3; iS++) {
                        if (!fs_in_delaunay_queue({iF, iS})) {
                            delaunay_check_queue.push({iF, iS});
                            in_delaunay_queue(iF, iS) = true;
                        }
                    }
                }
            }
            continue;
        }

        // If we've reached this point, the circumradius queue must be empty.
        // Now we double check to make sure we didn't miss anything
        // (which can happen rarely due to numerics), but don't do this more than a few times,
        // to avoid getting stuck in an infinite loop when numerical ultra-badness happens
        if (recheck_count < MAX_RECHECK_COUNT) {
            recheck_count++;
            bool any_found = false;
            if (delaunay_check_queue.empty() && circumradius_check_queue.empty()) {
                for (int iF = 0; iF < F.rows(); iF++) {
                    if (needs_refinement(iF)) {
                        circumradius_check_queue.push(make_tuple(area_weight(iF), face_area(iF), iF));
                        any_found = true;
                    }
                    for (int iS = 0; iS < 3; iS++) {
                        Vector2i fs{iF, iS};
                        if (!is_delaunay(G,l,fs) && !is_boundary_face_side(G,fs)
                            && !fs_in_delaunay_queue(fs)) {

                            delaunay_check_queue.push(fs);
                            in_delaunay_queue(iF, iS) = true;
                            any_found = true;
                        }
                    }
                }
            }

            if (!any_found) {
                // makes sure we don't recheck multiple times in a row
                break;
            }
        }
    } while (!delaunay_check_queue.empty() || !circumradius_check_queue.empty()
             || recheck_count < MAX_RECHECK_COUNT);

}

