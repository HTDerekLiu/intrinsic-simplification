#include "always_flip_to_ear.h"

static double face_side_score(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs)
{
    // compute the flip score for each face side
    using namespace global_variables;
    using namespace std;
    int vi, vj;
    get_face_side_vertices(F,fs,vi,vj);
    if (vi == vj){
        // SELF_EDGE_SCORE will be larger than opposite_corner_angle sum, so it will get flipped first
        return SELF_EDGE_SCORE; 
    }
    else
        return opposite_corner_angle(l, fs) + opposite_corner_angle(l, twin(G,fs)); 
}

template <typename T, typename A>
static int argmax(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

void always_flip_to_ear(
    const int & v,
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
    using namespace global_variables;

    assert(is_boundary_vertex(G,v2fs,v));
    assert (abs(gaussian_curvature_at_vertex(G,l,v2fs,v)) < 1e-7);
    assert(is_one_ring_fs_valid(F,G,v2fs,v));

    // get number of one-ring vertices and one ring fs
    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G,v2fs,v,fs_list);

    // break if this is an ear already
    if (fs_list.size() == 1)
        return;

    // compute the diamond angle sum for each one-ring face side
    fs_list.erase(fs_list.begin()); // remove the first face side because that is the boundary fs and it cannot be flipped

    vector<double> score_array;
    score_array.reserve(fs_list.size());
    for (int ii=0; ii<fs_list.size(); ii++)
    {
        Vector2i fs = fs_list[ii];
        score_array.push_back(face_side_score(F, G, l, fs));
    }

    // start flipping
    int idx, idx_next, idx_prev;
    while (true)
    {
        assert(is_one_ring_fs_valid(F,G,v2fs,v));

        // get a flippable face side
        idx = argmax(score_array);
        Vector2i fs_flip = fs_list[idx];

        // check whether fs_flip is a self edge
        bool is_self_edge =  (score_array[idx] == SELF_EDGE_SCORE);

        if (is_diamond_convex(G,l,fs_flip))
        {
            // flip the edge 
            map<pair<int, int>, pair<int, int>> fs_dict;
            flip_edge(fs_flip, F,G,l,A,v2fs,BC,F2V,fs_dict);

            // update fs list and ang array
            if ((idx+1) < fs_list.size()) // if has idx_next
            {
                idx_next = idx + 1;

                // relabel idx_next
                auto fs_next_before = make_pair(fs_list[idx_next](0), fs_list[idx_next](1));
                auto fs_next_after = fs_dict.find(fs_next_before)->second;
                fs_list[idx_next] << fs_next_after.first, fs_next_after.second;

                // update score array
                score_array[idx_next] = face_side_score(F, G, l, fs_list[idx_next]);
            }
            if (idx > 0) // if has idx_prev
            {
                idx_prev = idx - 1;

                // relabel idx_prev
                fs_list[idx_prev] = next(next(fs_list[idx]));

                // update score array
                score_array[idx_prev] = face_side_score(F, G, l, fs_list[idx_prev]);
            }

            // erase elements (TODO: can we do this without removing elements from vectors so that it is faster??)
            fs_list.erase(fs_list.begin()+idx);
            score_array.erase(score_array.begin()+idx);

            // delete twin of fs_flip if it is an edge connecting the same vertex
            if (is_self_edge){
                Vector2i fs_flip_twin = twin(G,fs_flip);
                for (int ii=0; ii<fs_list.size(); ii++){
                    if (is_same_face_side(fs_flip_twin, fs_list[ii])){

                        // update fs list and ang array
                        if ((ii+1) < fs_list.size()) // if has idx_next
                        {
                            int ii_next = ii + 1;

                            // relabel idx_next
                            auto fs_next_before = make_pair(fs_list[ii_next](0), fs_list[ii_next](1));
                            auto fs_next_after = fs_dict.find(fs_next_before)->second;
                            fs_list[ii_next] << fs_next_after.first, fs_next_after.second;

                            // update score array
                            score_array[ii_next] = face_side_score(F, G, l, fs_list[ii_next]);
                        }
                        if (ii > 0) // if has idx_prev
                        {
                            int ii_prev = ii - 1;

                            // relabel idx_prev
                            fs_list[ii_prev] = next(next(fs_list[ii]));

                            // update score array
                            score_array[ii_prev] = face_side_score(F, G, l, fs_list[ii_prev]);
                        }

                        // remove it
                        fs_list.erase(fs_list.begin()+ii);
                        score_array.erase(score_array.begin()+ii);
                        break;
                    }
                }
            }

            // stopping criteria
            if (fs_list.size() == 0)
                return;
        }
        else // this edge cannot be flipped due to non diamond
            score_array[idx] = 0.0;

        {
            // check whether ang_array is all zeros
            bool is_all_zeros = true;
            for (int ii=0; ii<score_array.size(); ii++)
                if (score_array[ii] > EPS)
                {
                    is_all_zeros = false;
                    break;
                }
            if (is_all_zeros)
                throw std::runtime_error("[Error] failed to flip to ear in always_flip_to_ear");
        }
    }
}
