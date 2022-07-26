#include "dynamical_system.hpp"

Eigen::VectorXd dynamical_system::func([[maybe_unused]] double t,
                                       const Eigen::VectorXd &x) {
  Eigen::VectorXd ret(xdim);

  ret(0) = (((1.0 / 2.0) * tanh(x[1]) + 1.0 / 2.0) * p[3] -
            p[2] * pow(x[0], 3) + x[0] - x[2]) *
           p[4];
  ret(1) = (-((1.0 / 2.0) * tanh(x[0]) + 1.0 / 2.0) * p[3] -
            p[2] * pow(x[1], 3) + x[1] - x[3]) *
           p[4];
  ret(2) = p[0] - p[1] * x[2] + x[0];
  ret(3) = p[0] - p[1] * x[3] + x[1];

  return ret;
}