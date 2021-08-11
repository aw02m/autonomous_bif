#include "ds_derivatives.hpp"

Eigen::VectorXd f(double t, const Eigen::VectorXd &x,
                  const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim);

  // rossler
  ret(0) = -x(1) - x(2);
  ret(1) = x(0) + ds.params(0) * x(1);
  ret(2) = ds.params(1) * x(0) - ds.params(2) * x(2) + x(0) * x(2);

  // // bvp
  // ret(0) = -x(2) + ds.params(0) * x(0) + std::tanh(ds.params(1) * x(0));
  // ret(1) = x(2) - ds.params(2) * x(1);
  // ret(2) = x(0) - x(1);

  return ret;
}

Eigen::MatrixXd dfdx(const Eigen::VectorXd &x, const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.xdim, ds.xdim);

  // rossler
  ret(0, 0) = 0.0;
  ret(0, 1) = -1.0;
  ret(0, 2) = -1.0;
  ret(1, 0) = 1.0;
  ret(1, 1) = ds.params(0);
  ret(1, 2) = 0.0;
  ret(2, 0) = ds.params(1) + x(2);
  ret(2, 1) = 0.0;
  ret(2, 2) = -ds.params(2) + x(0);

  // bvp
  // ret(0, 0) =
  //     ds.params(0) + 1/std::pow(std::cosh(ds.params(1) * x(0)), 2) *
  //     ds.params(1);
  // ret(0, 1) = 0;
  // ret(0, 2) = -1;
  // ret(1, 0) = 0;
  // ret(1, 1) = -ds.params(2);
  // ret(1, 2) = 1;
  // ret(2, 0) = 1;
  // ret(2, 1) = -1;
  // ret(2, 2) = 0;

  return ret;
}

Eigen::VectorXd dfdlambda(const Eigen::VectorXd &x,
                          const dynamical_system &ds) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(ds.xdim);

  switch (ds.var_param) {
  case 0:
    ret(1) = x(1);
    break;
  case 1:
    ret(2) = x(0);
    break;
  case 2:
    ret(2) = -x(2);
    break;
  }

  return ret;
}

std::vector<Eigen::MatrixXd> dfdxdx(const Eigen::VectorXd &x,
                                    const dynamical_system &ds) {
  std::vector<Eigen::MatrixXd> ret(ds.xdim,
                                   Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));

  // for (int i = 0; i < ds.xdim; i++) {
    // ret[i] = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
  // }

  ret[0](2, 2) = 1;
  ret[2](2, 0) = 1;
  // switch (k) {
  // case 0:
  //   ret(2, 2) = 1;
  //   break;
  // case 1:
  //   break;
  // case 2:
  //   ret(2, 0) = 1;
  //   break;
  // }

  return ret;
}

Eigen::MatrixXd dfdxdlambda(const Eigen::VectorXd &x,
                            const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);

  switch (ds.var_param) {
  case 0:
    ret(1, 1) = 1;
    break;
  case 1:
    ret(2, 0) = 1;
    break;
  case 2:
    ret(2, 2) = -1;
    break;
  }

  return ret;
}