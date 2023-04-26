#include "is_diamond_looped.h"

bool is_diamond_looped(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & s0)
{
    using namespace Eigen;

    Vector2i t0 = twin(G, s0);
	Vector2i t1 = next(t0);
	Vector2i t2 = next(t1);

	Vector2i s1 = next(s0);
	Vector2i s2 = next(s1);

	int	v0 = F(s0(0), s0(1));
	int	v1 = F(s1(0), s1(1));
	int	v2 = F(s2(0), s2(1));
	int	v3 = F(t2(0), t2(1));

    std::set<int> vs;
    vs.insert(v0);
    vs.insert(v1);
    vs.insert(v2);
    vs.insert(v3);
    return (vs.size() == 4);
}