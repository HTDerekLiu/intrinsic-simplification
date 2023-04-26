#include "remove_vector_element.h"

void remove_vector_element(
  const int idx,
  Eigen::VectorXd & vec)
  {
    unsigned int output_size = vec.size() - 1;

    if( idx < output_size )
        vec.segment(idx,output_size-idx) = vec.segment(idx+1,output_size-idx);

    vec.conservativeResize(output_size);
  }

void remove_vector_element(
  const int idx,
  Eigen::VectorXi & vec)
  {
    unsigned int output_size = vec.size() - 1;

    if( idx < output_size )
        vec.segment(idx,output_size-idx) = vec.segment(idx+1,output_size-idx);

    vec.conservativeResize(output_size);
  }
