#include "always_flip_to_degree_three.h"

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

static int argmax_vec(
    const Eigen::VectorXd & vec)
{
    // return the index of the max element in a vector
    double max_val = vec(0);
    int max_idx = 0;
    for (int ii=1; ii<vec.size(); ii++)
    {
        if (vec(ii) > max_val)
        {
            max_idx = ii;
            max_val = vec(ii);
        }
    }
    return max_idx;
}

static void relabel(
    const int v,
    const Eigen::MatrixXi & F,
    const int & idx,
    const std::map<std::pair<int, int>, std::pair<int, int>> & fs_dict,
    std::vector<Eigen::Vector2i> & fs_list)
{
    using namespace std;
    int degree = fs_list.size();

    // updateh idx_next
    {
        int idx_next = (idx+1) % degree;
        int f = fs_list[idx_next](0);
        int s = fs_list[idx_next](1);
        if (F(f, s) != v){ 
            // relabel the face side only if the starting vertex is not v
            // this will be the case for the neighboring edges of an self-edge 
            auto fs_next_before = make_pair(f, s);
            auto fs_next_after = fs_dict.find(fs_next_before)->second;
            fs_list[idx_next] << fs_next_after.first, fs_next_after.second;
        }
    }

    // update idx_prev
    {
        int idx_prev = (idx-1+degree) % degree;
        int f = fs_list[idx_prev](0);
        int s = fs_list[idx_prev](1);
        if (F(f, s) != v){ 
            // relabel the face side only if the starting vertex is not v
            // this will be the case for the neighboring edges of an self-edge 
            fs_list[idx_prev] = next(next(fs_list[idx]));
        }
    }
}

void always_flip_to_degree_three(
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

    assert (is_interior_vertex(G,v2fs,v));
    assert (abs(gaussian_curvature_at_vertex(G,l,v2fs,v)) < 1e-7);
    assert (is_one_ring_fs_valid(F,G,v2fs,v));

    // get number of one-ring vertices and one ring fs
    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G,v2fs,v,fs_list);
    int n = fs_list.size();

    // if already degree three, then don't do anything
    if (n == 3)
        return;

    // for each face side, we compute the score of a diamond
    VectorXd score_array(n);
    for (int ii=0; ii<n; ii++)
    {
        Vector2i fs = fs_list[ii];
        score_array(ii) = face_side_score(F,G,l,fs_list[ii]);
    }

    // start flipping to degree three
    while (1)
    {
        assert(is_one_ring_fs_valid(F,G,v2fs,v));

        int idx = argmax_vec(score_array);
        Vector2i fs_flip = fs_list[idx];

        // check whether fs_flip is a self edge
        bool is_self_edge =  (score_array(idx) == SELF_EDGE_SCORE);

        if (is_diamond_convex(G,l,fs_flip))
        {
            // flip the edge 
            map<pair<int, int>, pair<int, int>> fs_dict;
            flip_edge(fs_flip, F,G,l,A,v2fs,BC,F2V,fs_dict);

            // relabel the fs_list (to handle self edges)
            relabel(v, F, idx, fs_dict, fs_list);

            // update score array
            int idx_next = (idx+1) % fs_list.size();
            int idx_prev = (idx-1+fs_list.size()) % fs_list.size();
            score_array(idx_next) = face_side_score(F,G,l,fs_list[idx_next]);
            score_array(idx_prev) = face_side_score(F,G,l,fs_list[idx_prev]);

            // delete elements
            fs_list.erase(fs_list.begin()+idx);
            remove_vector_element(idx, score_array);

            // delete twin of fs_flip if it is an edge connecting the same vertex
            if (is_self_edge){
                Vector2i fs_flip_twin = twin(G,fs_flip);
                for (int ii=0; ii<fs_list.size(); ii++){
                    if (is_same_face_side(fs_flip_twin, fs_list[ii])){

                        // relabel the fs_list (to handle self edges)
                        relabel(v, F, ii, fs_dict, fs_list);

                        // update score array
                        int ii_next = (ii+1) % fs_list.size();
                        int ii_prev = (ii-1+fs_list.size()) % fs_list.size();
                        score_array(ii_next) = face_side_score(F,G,l,fs_list[ii_next]);
                        score_array(ii_prev) = face_side_score(F,G,l,fs_list[ii_prev]);

                        // remove it
                        fs_list.erase(fs_list.begin()+ii);
                        remove_vector_element(ii, score_array);
                        break;
                    }
                }
            }

            // stopping criteria
            if (fs_list.size() == 3){
                return;
            }
        }
        else // this edge cannot be flipped due to non diamond
            score_array(idx) = 0.0;
        
        if (score_array.isZero())
            throw std::runtime_error("[Error] failed to flip to degree three in always_flip_to_degree_three");
    }
    
}
