#include "signed_angle.h"

double signed_angle(
	const Eigen::Vector2d & u,
	const Eigen::Vector2d & v)
{
	using namespace std;
	double angle = atan2(v(1),v(0)) - atan2(u(1),u(0));
	if (angle < 0){
		angle += 2 * M_PI;
	}
	return angle;
}