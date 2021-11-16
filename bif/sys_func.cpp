#include "dynamical_system.hpp"

void dynamical_system::operator()(const Eigen::VectorXd &x,
                                Eigen::VectorXd &dxdt, const double /*t*/) {
  // rossler
  f(0) = -x(1) - x(2);
  f(1) = x(0) + p(0) * x(1);
  f(2) = p(1) * x(0) - p(2) * x(2) + x(0) * x(2);

  dfdx(0, 0) = 0.0;
  dfdx(0, 1) = -1.0;
  dfdx(0, 2) = -1.0;
  dfdx(1, 0) = 1.0;
  dfdx(1, 1) = p(0);
  dfdx(1, 2) = 0.0;
  dfdx(2, 0) = p(1) + x(2);
  dfdx(2, 1) = 0.0;
  dfdx(2, 2) = -p(2) + x(0);

  if (mode != 0) {
    dfdlambda = Eigen::VectorXd::Zero(xdim);
    switch (var_param) {
    case 0:
      dfdlambda(1) = x(1);
      break;
    case 1:
      dfdlambda(2) = x(0);
      break;
    case 2:
      dfdlambda(2) = -x(2);
      break;
    }

    dfdxdx =
        std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));
    dfdxdx[0](2, 2) = 1;
    dfdxdx[2](2, 0) = 1;

    dfdxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);
    switch (var_param) {
    case 0:
      dfdxdlambda(1, 1) = 1;
      break;
    case 1:
      dfdxdlambda(2, 0) = 1;
      break;
    case 2:
      dfdxdlambda(2, 2) = -1;
      break;
    }
  }

  /******************************************************/
  /*     DO NOT EDIT BELOW (variational equations)      */
  /******************************************************/
  unsigned int counter = xdim;
  // variational state (transform to matrix shape for easy producting)
  Eigen::MatrixXd state_dphidx = x(Eigen::seqN(counter, size_dphidx));
  state_dphidx.resize(xdim, xdim);
  counter += size_dphidx;
  Eigen::VectorXd state_dphidlambda = x(Eigen::seqN(counter, size_dphidlambda));
  counter += size_dphidlambda;
  std::vector<Eigen::MatrixXd> state_dphidxdx(
      xdim, Eigen::MatrixXd::Zero(xdim, xdim));
  Eigen::MatrixXd temp;
  for (int i = 0; i < xdim; i++) {
    temp = x(Eigen::seqN(counter + size_dphidx * i, size_dphidx));
    temp.resize(xdim, xdim);
    state_dphidxdx[i] = temp;
    temp.resize(size_dphidx, 1);
  }
  counter += size_dphidxdx;
  Eigen::MatrixXd state_dphidxdlambda =
      x(Eigen::seqN(counter, size_dphidxdlambda));

  counter = 0;

  // phi
  dxdt(Eigen::seqN(counter, xdim)) = f;
  counter += xdim;

  // dphidx
  dxdt(Eigen::seqN(counter, size_dphidx)) = dfdx * state_dphidx;
  counter += size_dphidx;

  if (mode != 0) {
    // dphidlambda
    dxdt(Eigen::seqN(counter, xdim)) = dfdx * state_dphidlambda + dfdlambda;
    counter += size_dphidlambda;

    // dphidxdx
    for (int i = 0; i < xdim; i++) {
      dxdt(Eigen::seqN(counter, size_dphidx)) =
          dfdxdx[i] * state_dphidx * state_dphidx + dfdx * state_dphidxdx[i];
      counter += size_dphidx;
    }

    // dphidxdlambda
    for (int i = 0; i < xdim; i++) {
      dxdt(Eigen::seqN(counter, xdim)) =
          dfdxdx[i] * state_dphidlambda * state_dphidx +
          dfdxdlambda * state_dphidx.col(i) + dfdx * state_dphidxdlambda.col(i);
      counter += xdim;
    }
  }
}