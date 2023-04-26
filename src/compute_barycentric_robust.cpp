#include "compute_barycentric_robust.h"

Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c)
{
    bool is_degenerated, is_inside;
    double distance_to_valid;
    return compute_barycentric_robust(p, a, b, c, is_degenerated, is_inside, distance_to_valid);
}

Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c,
    const double & degenerate_threshold)
{
    bool is_degenerated, is_inside;
    double distance_to_valid;
    return compute_barycentric_robust(p, a, b, c, degenerate_threshold, is_degenerated, is_inside, distance_to_valid);
}

Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c,
    bool & is_degenerated,
    bool & is_inside,
    double & distance_to_valid)
{
    double degenerate_threshold = 1e-5;
    return compute_barycentric_robust(p, a, b, c, degenerate_threshold, is_degenerated, is_inside, distance_to_valid);
}


Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c,
    const double & degenerate_threshold,
    bool & is_degenerated,
    bool & is_inside,
    double & distance_to_valid)
{
    using namespace Eigen;
    using namespace std;
    using namespace global_variables;

    Vector3d bc;

    // build system matrix
    MatrixXd A(3, 4);
    A << a(0), b(0), c(0), p(0),
         a(1), b(1), c(1), p(1),
            1,    1,    1,    1;

    // compute barycentric coordinates in the barycentric coordinates
    Vector4d raw_b;
    {
        Matrix3d tmp;
        
        tmp << A(0,1), A(0,2), A(0,3),
               A(1,1), A(1,2), A(1,3),
               A(2,1), A(2,2), A(2,3);
        raw_b(0) = tmp.determinant();

        tmp << A(0,0), A(0,2), A(0,3),
               A(1,0), A(1,2), A(1,3),
               A(2,0), A(2,2), A(2,3);
        raw_b(1) = -tmp.determinant();

        tmp << A(0,0), A(0,1), A(0,3),
               A(1,0), A(1,1), A(1,3),
               A(2,0), A(2,1), A(2,3);
        raw_b(2) = tmp.determinant();
        
        tmp << A(0,0), A(0,1), A(0,2),
               A(1,0), A(1,1), A(1,2),
               A(2,0), A(2,1), A(2,2);
        raw_b(3) = -tmp.determinant();
    }

    // determine whether it is degenerated
    double epsilon = degenerate_threshold;
    is_degenerated = abs(raw_b(3)) < epsilon;

    // determind whether this bary is inside
    {
        if (raw_b(3) > 0)
        {
            is_inside = (
                -epsilon <= -raw_b(0) &&
                -epsilon <= -raw_b(1) &&
                -epsilon <= -raw_b(2) &&
                -raw_b(0) <= (raw_b(3)+epsilon) &&
                -raw_b(1) <= (raw_b(3)+epsilon) &&
                -raw_b(2) <= (raw_b(3)+epsilon) 
            );
        }
        else
        {
            is_inside = (
                (raw_b(3)-epsilon) <= -raw_b(0) &&
                (raw_b(3)-epsilon) <= -raw_b(1) &&
                (raw_b(3)-epsilon) <= -raw_b(2) &&
                -raw_b(0) <= epsilon &&
                -raw_b(1) <= epsilon &&
                -raw_b(2) <= epsilon 
            );
        }
    }

    // compute output barycentric coordinates
    if (!is_degenerated && is_inside) // valid case
    {
        // get barycentric
        bc << -raw_b.head(3);
        bc = bc.array() / raw_b(3);

        // compute distance to valid barycentric
        {
            Vector3d dist2zero, dist2one;
            dist2zero = -bc;
            for (int ii=0; ii<dist2zero.size(); ii++)
                if (dist2zero(ii) < 0)
                    dist2zero(ii) = 0.0;
            dist2one = bc.array() - 1.0;
            for (int ii=0; ii<dist2one.size(); ii++)
                if (dist2one(ii) < 0)
                    dist2one(ii) = 0.0;

            distance_to_valid = max(dist2zero.maxCoeff(), dist2one.maxCoeff());
        }

        // clip to valid barycentric (some of them are caused by numerical issues)
        for (int ii = 0; ii < bc.size(); ii++)
        {
            if (bc(ii) < 0)
                bc(ii) = 0;
            else if (bc(ii) > 1)
                bc(ii) = 1;
        }

        // normalize to have sum 1
        bc = bc.array() / bc.sum();
    }
    else if (!is_degenerated && !is_inside) // outside triangle
    {
        // get barycentric
        bc << -raw_b.head(3);
        bc = bc.array() / raw_b(3);

        // compute distance to valid barycentric
        {
            Vector3d dist2zero, dist2one;
            dist2zero = -bc;
            for (int ii=0; ii<dist2zero.size(); ii++)
                if (dist2zero(ii) < 0)
                    dist2zero(ii) = 0.0;
            dist2one = bc.array() - 1.0;
            for (int ii=0; ii<dist2one.size(); ii++)
                if (dist2one(ii) < 0)
                    dist2one(ii) = 0.0;

            distance_to_valid = max(dist2zero.maxCoeff(), dist2one.maxCoeff());
        }
    }
    else // degenerated triangle
    {
        bc << DOUBLE_INF, DOUBLE_INF, DOUBLE_INF;
        distance_to_valid = DOUBLE_INF;
    }

    return bc;
}