#ifndef SIGNED_ANGLE
#define SIGNED_ANGLE

#include <Eigen/Core>
#include <math.h>

#include <pi.h>

// Given two 2d vectors (u, v), compute the signed angle between them
double signed_angle(
	const Eigen::Vector2d & u,
	const Eigen::Vector2d & v);
#endif
