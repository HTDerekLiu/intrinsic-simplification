#ifndef TWIN
#define TWIN

#include <Eigen/Core>
#include <Eigen/Dense>

// For a given face-side fs, returns the twin face-side 
// Inputs
//     G: |F|x6 gluing map G,
//     fs: a face-side tuple (f,s)
// Outputs
//     fs_twin The twin of face-side fs
Eigen::Vector2i twin(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs);

#endif