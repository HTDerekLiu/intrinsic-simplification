#include "bake_texture.h"

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
    MatrixXi colors(24, 3); // from colorbrewer2.org
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


    vector<pair<int, int>> coord_faces;
	for (int fIdx=0; fIdx<F2V.size(); fIdx++){
		if (F2V[fIdx].size() > 0){
			for (auto idx : F2V[fIdx]){
                if (idx >= nV)
                    coord_faces.emplace_back(idx - nV, fIdx);
			}
		}
	}
    std::sort(coord_faces.begin(), coord_faces.end());

    int tex_size = hit_mask.size();
    int tex_width = (int) sqrt(tex_size);
    VectorXi query_point_indices(tex_size);
    query_point_indices.setConstant(-1);
    int coord_idx = 0;
    for (int i = 0; i < tex_size; i++) {
        int idx = i + nV;
        if (hit_mask(i)) {
            query_point_indices(i) = coord_faces[coord_idx++].second;
        }
    }

    VectorXi query_point_indices_copy(query_point_indices);
    // fill in blank cells with an adjacent color if possible
    for (int i = 0; i < tex_width; i++) {
        for (int j = 0; j < tex_width; j++) {
        int coord = i * tex_width + j;
        if (query_point_indices(coord) != -1) continue;
        int xoff, yoff;
        xoff = -1; yoff = 0;
        if (j+xoff >= 0 && j+xoff < tex_width && i+yoff >= 0 && i+yoff < tex_width && query_point_indices((i+yoff)*tex_width+j+xoff) != -1) {
            query_point_indices_copy(coord) = query_point_indices((i+yoff)*tex_width+j+xoff);
            continue;
        }
        xoff = 1;
        if (j+xoff >= 0 && j+xoff < tex_width && i+yoff >= 0 && i+yoff < tex_width && query_point_indices((i+yoff)*tex_width+j+xoff) != -1) {
            query_point_indices_copy(coord) = query_point_indices((i+yoff)*tex_width+j+xoff);
            continue;
        }
        xoff = 0; yoff = 1;
        if (j+xoff >= 0 && j+xoff < tex_width && i+yoff >= 0 && i+yoff < tex_width && query_point_indices((i+yoff)*tex_width+j+xoff) != -1) {
            query_point_indices_copy(coord) = query_point_indices((i+yoff)*tex_width+j+xoff);
            continue;
        }
        yoff = -1;
        if (j+xoff >= 0 && j+xoff < tex_width && i+yoff >= 0 && i+yoff < tex_width && query_point_indices((i+yoff)*tex_width+j+xoff) != -1) {
            query_point_indices_copy(coord) = query_point_indices((i+yoff)*tex_width+j+xoff);
            continue;
        }
        }
    }
    query_point_indices = query_point_indices_copy;

    int maxmax = 0;
    int nF = F2V.size();
    VectorXi query_point_color_indices(nF);
    query_point_color_indices.setConstant(-1);
    for (int f = 0; f < nF; f++) {
        // if (f%100 == 0)
        //     cout << "query_point_color_indices: " << f << "/" << nF << endl;

        set<int> f1_vs = {F(f, 0), F(f, 1), F(f, 2)};
        set<int> nbrs;
        if (*f1_vs.begin() < 0) continue;
        for (int f2 = 0; f2 < nF; f2++) {
            if (f == f2) continue;
            int cnt = f1_vs.count(F(f2, 0)) + f1_vs.count(F(f2, 1)) + f1_vs.count(F(f2, 2));
            if (cnt >= 1) {
                if (query_point_color_indices(f2) != -1)
                    nbrs.insert(query_point_color_indices(f2));
            }
        }

        int mincoord = 0;
        while (nbrs.count(mincoord) > 0) mincoord++;

        if (mincoord >= colors.rows())
            throw std::runtime_error("[Error] not enough color in back texture cpp");

        query_point_color_indices(f) = mincoord;
    }


    MatrixXi query_point_colors(tex_size, 3);
    for (int i = 0; i < tex_size; i++) {
        if (query_point_indices(i) == -1 || query_point_color_indices(query_point_indices(i)) == -1)
            query_point_colors.row(i) = Vector3i{0, 0, 0};
        else
            query_point_colors.row(i) = colors.row(query_point_color_indices(query_point_indices(i)));
    }

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
