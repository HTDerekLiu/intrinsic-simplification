#include "ray_line_intersection.h"
 
void ray_line_intersection(
	const Eigen::Vector2d & o,
	const Eigen::Vector2d & d,
	const Eigen::Vector2d & v1,
	const Eigen::Vector2d & v2,
	double & length,
	double & ratio,
	double & angle)
{
    using namespace std;
	using namespace Eigen;
	using namespace global_variables;

	Vector2d u1 = o - v1;
	Vector2d u2 = v2 - v1;
	Vector2d u3 = {-d(1), d(0)};
	double u2_cross_u1 = u2(0)*u1(1) - u1(0)*u2(1);
	length = u2_cross_u1 / u2.dot(u3); // length along d
	ratio = u1.dot(u3)/ u2.dot(u3); // intersection is: v1 + ratio*(v2-v1)
	angle = signed_angle(d, v2-v1); // angle is always on the RHS of the path

	if ((ratio < -1e-16) || (ratio > (1+1e-16))){
		length = DOUBLE_NAN;
		ratio = DOUBLE_NAN;
		angle = DOUBLE_NAN;
	}
	else{
		ratio = max(0.0, min(ratio, 1.0));
	}
}