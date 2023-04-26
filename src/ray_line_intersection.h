#ifndef RAY_LINE_INTERSECTION
#define RAY_LINE_INTERSECTION

#include <Eigen/Core>
#include <math.h>
#include "signed_angle.h"
#include "global_variables.h"

/*
    computes the ray "o+length*d" and line segment (v1, v2) interesection
    
    Inputs:
    o: (2,) array of ray origin
    d: (2,) array of ray direction
    v1: (2,) array of a end point of a line
    v2: (2,) array of the other end point of a line

    Outputs
    length: length from the origin o to the hit location
    ratio: ratio of the hit location so that p = v1 + ratio*(v2-v1)
    angle: angle between the direction and the line segment (see figure)

    Notes
    #         o 
    #         |    
    #         |   
    #         d    
    #         |     
    #         |      
    #     ang V       
    #  v1 ---------- v2
*/
void ray_line_intersection(
	const Eigen::Vector2d & o,
	const Eigen::Vector2d & d,
	const Eigen::Vector2d & v1,
	const Eigen::Vector2d & v2,
	double & length,
	double & ratio,
	double & angle);
#endif