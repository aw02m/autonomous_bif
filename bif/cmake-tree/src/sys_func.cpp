#include "dynamical_system.hpp"

void dynamical_system::sys_func(const Eigen::VectorXd &x, const double /*t*/) {
  f(0) = -x[1] - x[2];
  f(1) = p[0] * x[1] + x[0];
  f(2) = p[1] * x[0] - p[2] * x[2] + x[0] * x[2];

  dfdx(0, 0) = 0;
  dfdx(0, 1) = -1;
  dfdx(0, 2) = -1;
  dfdx(1, 0) = 1;
  dfdx(1, 1) = p[0];
  dfdx(1, 2) = 0;
  dfdx(2, 0) = p[1] + x[2];
  dfdx(2, 1) = 0;
  dfdx(2, 2) = -p[2] + x[0];

  switch (var_param) {
  case 0:
    dfdlambda(0) = 0;
    dfdlambda(1) = x[1];
    dfdlambda(2) = 0;
    break;
  case 1:
    dfdlambda(0) = 0;
    dfdlambda(1) = 0;
    dfdlambda(2) = x[0];
    break;
  case 2:
    dfdlambda(0) = 0;
    dfdlambda(1) = 0;
    dfdlambda(2) = -x[2];
    break;
  }

  dfdxdx[0](0, 0) = 0;
  dfdxdx[0](0, 1) = 0;
  dfdxdx[0](0, 2) = 0;
  dfdxdx[0](1, 0) = 0;
  dfdxdx[0](1, 1) = 0;
  dfdxdx[0](1, 2) = 0;
  dfdxdx[0](2, 0) = 0;
  dfdxdx[0](2, 1) = 0;
  dfdxdx[0](2, 2) = 1;
  dfdxdx[1](0, 0) = 0;
  dfdxdx[1](0, 1) = 0;
  dfdxdx[1](0, 2) = 0;
  dfdxdx[1](1, 0) = 0;
  dfdxdx[1](1, 1) = 0;
  dfdxdx[1](1, 2) = 0;
  dfdxdx[1](2, 0) = 0;
  dfdxdx[1](2, 1) = 0;
  dfdxdx[1](2, 2) = 0;
  dfdxdx[2](0, 0) = 0;
  dfdxdx[2](0, 1) = 0;
  dfdxdx[2](0, 2) = 0;
  dfdxdx[2](1, 0) = 0;
  dfdxdx[2](1, 1) = 0;
  dfdxdx[2](1, 2) = 0;
  dfdxdx[2](2, 0) = 1;
  dfdxdx[2](2, 1) = 0;
  dfdxdx[2](2, 2) = 0;

  switch (var_param) {
  case 0:
    dfdxdlambda(0, 0) = 0;
    dfdxdlambda(0, 1) = 0;
    dfdxdlambda(0, 2) = 0;
    dfdxdlambda(1, 0) = 0;
    dfdxdlambda(1, 1) = 1;
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(2, 0) = 0;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = 0;
    break;
  case 1:
    dfdxdlambda(0, 0) = 0;
    dfdxdlambda(0, 1) = 0;
    dfdxdlambda(0, 2) = 0;
    dfdxdlambda(1, 0) = 0;
    dfdxdlambda(1, 1) = 0;
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(2, 0) = 1;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = 0;
    break;
  case 2:
    dfdxdlambda(0, 0) = 0;
    dfdxdlambda(0, 1) = 0;
    dfdxdlambda(0, 2) = 0;
    dfdxdlambda(1, 0) = 0;
    dfdxdlambda(1, 1) = 0;
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(2, 0) = 0;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = -1;
    break;
  }
}