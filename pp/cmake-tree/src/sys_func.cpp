#include "dynamical_system.hpp"

Eigen::VectorXd dynamical_system::func([[maybe_unused]] double t,
                                       const Eigen::VectorXd &x) {
  Eigen::VectorXd ret(xdim);

  ret(0) = -x[1] - x[2];
  ret(1) = p[0] * x[1] + x[0];
  ret(2) = p[1] * x[0] - p[2] * x[2] + x[0] * x[2];

  return ret;
}