#include "dynamical_system.hpp"

Eigen::VectorXd dynamical_system::func([[maybe_unused]] double t,
                                       const Eigen::VectorXd &x) {
  Eigen::VectorXd ret(xdim);

  ret(0) = tanh(p[1] * x[0]) + p[0] * x[0] - x[2];
  ret(1) = x[2] - x[1] / p[2];
  ret(2) = x[0] - x[1];

  return ret;
}