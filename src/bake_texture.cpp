#include "bake_texture.h"

#define N_COLORS 24

void write_png(
    const std::string & png_path,
    unsigned char* image_buffer,
    int w,
    int h,
    int channels = 4)
{
    stbi_write_png(png_path.c_str(), w, h, channels, image_buffer, channels * w);
}

void bake_texture(
                  const std::string & png_path,
                  std::vector<unsigned char> & image_rgba_buffer) {
    int tex_size = image_rgba_buffer.size() / 4;
    int tex_width = (int) sqrt(tex_size);
    write_png(png_path, image_rgba_buffer.data(), tex_width, tex_width);
}

void bake_texture(
                  const std::string & png_path,
                  const Eigen::MatrixXi & F,
                  const std::vector<std::vector<int>> & F2V,
                  const Eigen::Matrix<bool, Eigen::Dynamic, 1> & hit_mask,
                  const int & nV)
{
    std::vector<unsigned char> image_rgba_buffer;
    bake_texture(image_rgba_buffer, F, F2V, hit_mask, nV);
    bake_texture(png_path, image_rgba_buffer);
}

void bake_texture(
    std::vector<unsigned char> & image_buffer,
    const Eigen::MatrixXi & F,
    const std::vector<std::vector<int>> & F2V,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> & hit_mask,
    const int & nV)
{
    using namespace std;
    using namespace Eigen;

    // face colors
    MatrixXi colors(N_COLORS, 3); // from colorbrewer2.org
    colors << 251, 128, 114,
              128, 177, 211,
              253, 180, 98,
              141, 211, 199,
              255, 255, 179,
              190, 186, 218,
              179, 222, 105,
              252, 205, 229,
              217, 217, 217,
              188, 128, 189,
              204, 235, 197,
              255, 237, 111,
              166,206, 227,
              31,120,180,
              178,223,138,
              51,160,44,
              251,154,153,
              227,26,28,
              253,191,111,
              255,127,0,
              202,178,214,
              106,61,154,
              255,255,153,
              177,89,40;


    //== identify which face each pixel lies in
    vector<pair<int, int>> coord_faces;
    for (int fIdx = 0; fIdx < F2V.size(); fIdx++){
        if (F2V[fIdx].size() > 0){
            for (auto idx : F2V[fIdx]){
                // idx - nV is the pixel index within the set of pixels where hit_mask is true
                if (idx >= nV) {
                    // store (pixel index, face index)
                    coord_faces.emplace_back(idx - nV, fIdx);
                }
            }
        }
    }
    // sort entries to lie in pixel order
    std::sort(coord_faces.begin(), coord_faces.end());

    int tex_size = hit_mask.size();
    int tex_width = (int) sqrt(tex_size);
    VectorXi query_point_indices(tex_size); // face index that pixels lie in
    query_point_indices.setConstant(-1);    // initialize to -1, pixels outside faces will remain -1
    int coord_idx = 0;
    for (int i = 0; i < tex_size; i++) {
        if (hit_mask(i)) {
            query_point_indices(i) = coord_faces[coord_idx++].second;
        }
    }

    //== color faces so that no two faces which share a vertex are the same color
    int n_coarse_vertices = F.maxCoeff() + 1; // note that nV is the number of _original_ vertices
    // find a list of which faces are incident on each vertex
    int nF = F2V.size();
    vector<vector<int>> vertex_faces(n_coarse_vertices, vector<int>{});
    for (int iF = 0; iF < nF; iF++) {
        for (int iC = 0; iC < 3; iC++) {
            if (F(iF, iC) < 0) continue;
            vertex_faces[F(iF, iC)].push_back(iF);
        }
    }
    VectorXi query_point_color_indices(nF);
    query_point_color_indices.setConstant(-1);
    for (int iF = 0; iF < nF; iF++) {
        array<bool, N_COLORS> used_color;
        used_color.fill(false);
        // iterate over neighbors of face f
        for (int iC = 0; iC < 3; iC++) {
            if (F(iF, iC) < 0) continue;
            for (int iN : vertex_faces[F(iF, iC)]) {
                if (iF == iN) continue;
                // if iN has been assigned a color, record it to ensure
                // that iF will pick a different color
                if (query_point_color_indices(iN) != -1)
                    used_color[query_point_color_indices(iN)] = true;
            }
        }

        // take the smallest color not already claimed by a neighbor of f
        int free_color = 0;
        while (used_color[free_color]) free_color++;

        if (free_color >= N_COLORS)
            throw std::runtime_error("[Error] not enough color in back texture cpp");

        query_point_color_indices(iF) = free_color;
    }

    //== record final pixel colors
    MatrixXi query_point_colors(tex_size, 3);
    for (int i = 0; i < tex_size; i++) {
        if (query_point_indices(i) == -1 || query_point_color_indices(query_point_indices(i)) == -1)
            query_point_colors.row(i) = Vector3i{0, 0, 0};
        else
            query_point_colors.row(i) = colors.row(query_point_color_indices(query_point_indices(i)));
    }

    // write colors to buffer
    image_buffer.resize(4 * tex_size);
    for (int j = 0; j < tex_width; j++) {
        for (int i = 0; i < tex_width; i++) {
            image_buffer[4 * (tex_width * j + i) + 0] = (unsigned char)(query_point_colors(j*tex_width+i, 0));
            image_buffer[4 * (tex_width * j + i) + 1] = (unsigned char)(query_point_colors(j*tex_width+i, 1));
            image_buffer[4 * (tex_width * j + i) + 2] = (unsigned char)(query_point_colors(j*tex_width+i, 2));
            image_buffer[4 * (tex_width * j + i) + 3] = (unsigned char)255;
        }
    }
}
