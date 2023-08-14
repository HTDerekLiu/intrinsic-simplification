#include "connected_components.h"

void connected_components(
                          const Eigen::MatrixXi & F,
                          const Eigen::MatrixXi & G,
                          int & n_components,
                          Eigen::VectorXi & v_id,
                          Eigen::VectorXi & f_id)
{
    using namespace Eigen;
    using namespace std;

    //== assign face components
    int iC = 0;
    f_id.resize(F.rows());
    f_id.setConstant(-1); // id of -1 indicates face is unvisited

    for (int iF = 0; iF < F.rows(); iF++) {
        if (f_id(iF) >= 0) continue;

        f_id(iF) = iC;

        queue<int> to_visit;
        to_visit.push(iF);

        while (!to_visit.empty()) {
            int iG = to_visit.front();
            to_visit.pop();

            // use glue map to iterate over neighbors
            for (size_t iS = 0; iS < 3; iS++) {
                int iH = G(iG, 2 * iS);

                // check if face has neighbor
                if (iH == global_variables::GHOST_INDEX) continue;

                // id of -1 indicates face is unvisited
                if (f_id(iH) < 0) {
                    f_id(iH) = iC;
                    to_visit.push(iH);
                }
            }
        }

        iC++;
    }

    //== assign vertex components
    int nV = F.maxCoeff() + 1;
    v_id.resize(nV);
    v_id.setConstant(-1); // id of -1 indicates vertex is unvisited

    // first, assign components to vertices contained in faces
    for (int iF = 0; iF < F.rows(); iF++) {
        for (int iS = 0; iS < 3; iS++) {
            v_id(F(iF, iS)) = f_id(iF);
        }
    }

    // then, assign components to any dangling vertices
    for (int iV = 0; iV < nV; iV++) {
        if (v_id(iV) < 0) v_id(iV) = iC++;
    }

    // record total number of assigned components
    n_components = iC;
}
