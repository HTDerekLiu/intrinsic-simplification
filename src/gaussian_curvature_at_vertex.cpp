#include "gaussian_curvature_at_vertex.h"

double gaussian_curvature_at_vertex(
    const Eigen::MatrixXi &G, 
    const Eigen::MatrixXd &l, 
    const Eigen::MatrixXi &v2fs, 
    const int v) 
{
  using namespace std;
  using namespace Eigen;
  using namespace global_variables;

  vector<Vector2i> fs_list;
  vertex_one_ring_face_sides(G, v2fs, v, fs_list);
  double angle_sum = 0.0;
  
  for (int ii=0; ii<fs_list.size(); ii++)
  {
    Vector2i fs = fs_list[ii];
    double angle = opposite_corner_angle(l, next(fs));
    if (isnan(angle))
    { 
      return DOUBLE_NAN;
    }
    angle_sum += angle;
  }

  if (is_interior_vertex(G,v2fs,v))
    return 2 * M_PI - angle_sum;
  else
    return M_PI - angle_sum;
}

