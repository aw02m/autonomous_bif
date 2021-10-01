#include "ds_derivatives.hpp"
#include "dynamical_system.hpp"

Eigen::VectorXd f(double t, const Eigen::VectorXd &x,
                  const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim);

  // // rossler
  // ret(0) = -x(1) - x(2);
  // ret(1) = x(0) + ds.params(0) * x(1);
  // ret(2) = ds.params(1) * x(0) - ds.params(2) * x(2) + x(0) * x(2);

  // // bvp
  ret(0) = -x(2) + ds.params(0) * x(0) + std::tanh(ds.params(1) * x(0));
  ret(1) = x(2) - ds.params(2) * x(1);
  ret(2) = x(0) - x(1);

  // sprott
  // ret(0) = x(1);
  // ret(1) = x(2) + ds.params(1);
  // ret(2) = -x(1) + 0.1 * x(0) * x(0) + 1.1 * x(0) * x(2) + ds.params(0);

  return ret;
}

Eigen::MatrixXd dfdx(const Eigen::VectorXd &x, const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.xdim, ds.xdim);

  // // rossler
  // ret(0, 0) = 0.0;
  // ret(0, 1) = -1.0;
  // ret(0, 2) = -1.0;
  // ret(1, 0) = 1.0;
  // ret(1, 1) = ds.params(0);
  // ret(1, 2) = 0.0;
  // ret(2, 0) = ds.params(1) + x(2);
  // ret(2, 1) = 0.0;
  // ret(2, 2) = -ds.params(2) + x(0);

  // bvp
  ret(0, 0) = ds.params(0) +
              1 / std::pow(std::cosh(ds.params(1) * x(0)), 2) * ds.params(1);
  ret(0, 1) = 0;
  ret(0, 2) = -1;
  ret(1, 0) = 0;
  ret(1, 1) = -ds.params(2);
  ret(1, 2) = 1;
  ret(2, 0) = 1;
  ret(2, 1) = -1;
  ret(2, 2) = 0;

  // sprott
  // ret(0, 0) = 0.0;
  // ret(0, 1) = 1.0;
  // ret(0, 2) = 0.0;
  // ret(1, 0) = 0.0;
  // ret(1, 1) = 0.0;
  // ret(1, 2) = 1.0;
  // ret(2, 0) = 0.2 * x(0) + 1.1 * x(2);
  // ret(2, 1) = -1.0;
  // ret(2, 2) = 1.1 * x(0);

  return ret;
}

Eigen::VectorXd dfdlambda(const Eigen::VectorXd &x,
                          const dynamical_system &ds) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(ds.xdim);

  // rossler
  // switch (ds.var_param) {
  // case 0:
  //   ret(1) = x(1);
  //   break;
  // case 1:
  //   ret(2) = x(0);
  //   break;
  // case 2:
  //   ret(2) = -x(2);
  //   break;
  // }

  // bvp
  switch (ds.var_param) {
  case 0:
    ret(0) = x(0);
    break;
  case 1:
    ret(0) = x(0) / std::pow(std::cosh(ds.params(1) * x(0)), 2);
    break;
  case 2:
    ret(1) = -x(1);
  }

  // sprott
  // switch (ds.var_param) {
  // case 0:
  //   ret(2) = 1;
  //   break;
  // case 1:
  //   ret(1) = 1;
  // }

  return ret;
}

std::vector<Eigen::MatrixXd> dfdxdx(const Eigen::VectorXd &x,
                                    const dynamical_system &ds) {
  std::vector<Eigen::MatrixXd> ret(ds.xdim,
                                   Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));

  // // rossler
  // ret[0](2, 2) = 1;
  // ret[2](2, 0) = 1;

  // bvp
  double b = ds.params(1);
  ret[0](0, 0) =
      (-2 * b * b * std::sinh(b * x(0))) / std::pow(std::cosh(b * x(0)), 3);

  // sprott
  // ret[0](2, 0) = 0.2;
  // ret[0](2, 2) = 1.1;
  // ret[2](2, 0) = 1.1;

  return ret;
}

Eigen::MatrixXd dfdxdlambda(const Eigen::VectorXd &x,
                            const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);

  // // rossler
  // switch (ds.var_param) {
  // case 0:
  //   ret(1, 1) = 1;
  //   break;
  // case 1:
  //   ret(2, 0) = 1;
  //   break;
  // case 2:
  //   ret(2, 2) = -1;
  //   break;
  // }

  // bvp
  double b = ds.params(1);
  switch (ds.var_param) {
  case 0:
    ret(0, 0) = 1;
    break;
  case 1:
    // ret(0, 0) = 1 / std::pow(std::cosh(b * x(0)), 2) +
    //             (-2 * b * x(0) * std::sinh(b * x(0))) /
    //                 (std::pow(std::cosh(b * x(0)), 3));
    // ret(0, 0) = (std::pow(std::cosh(b * x(0)), 2) -
    //              b * x(0) * std::sinh(2 * b * x(0))) /
    //             std::pow(std::cosh(b * x(0)), 4);
    ret(0, 0) = (std::cosh(b * x(0)) - 2 * b * x(0) * std::sinh(b * x(0))) /
                std::pow(std::cosh(b * x(0)), 3);
    break;
    case 2:
    ret(1, 1) = -1;
  }

  // sprott -> nothing (all elements are zero)

  return ret;
}