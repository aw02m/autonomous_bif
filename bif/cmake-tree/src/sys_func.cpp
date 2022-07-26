#include "dynamical_system.hpp"

void dynamical_system::operator()(const Eigen::VectorXd &x,
                                  Eigen::VectorXd &dxdt, const double /*t*/) {

  f(0) = (((1.0 / 2.0) * tanh(x[1]) + 1.0 / 2.0) * p[3] - p[2] * pow(x[0], 3) +
          x[0] - x[2]) *
         p[4];
  f(1) = (-((1.0 / 2.0) * tanh(x[0]) + 1.0 / 2.0) * p[3] - p[2] * pow(x[1], 3) +
          x[1] - x[3]) *
         p[4];
  f(2) = p[0] - p[1] * x[2] + x[0];
  f(3) = p[0] - p[1] * x[3] + x[1];

  dfdx(0, 0) = (-3 * p[2] * pow(x[0], 2) + 1) * p[4];
  dfdx(0, 1) = (1.0 / 2.0 - 1.0 / 2.0 * pow(tanh(x[1]), 2)) * p[3] * p[4];
  dfdx(0, 2) = -p[4];
  dfdx(0, 3) = 0;
  dfdx(1, 0) = ((1.0 / 2.0) * pow(tanh(x[0]), 2) - 1.0 / 2.0) * p[3] * p[4];
  dfdx(1, 1) = (-3 * p[2] * pow(x[1], 2) + 1) * p[4];
  dfdx(1, 2) = 0;
  dfdx(1, 3) = -p[4];
  dfdx(2, 0) = 1;
  dfdx(2, 1) = 0;
  dfdx(2, 2) = -p[1];
  dfdx(2, 3) = 0;
  dfdx(3, 0) = 0;
  dfdx(3, 1) = 1;
  dfdx(3, 2) = 0;
  dfdx(3, 3) = -p[1];

  switch (var_param) {
  case 0:
    dfdlambda(0) = 0;
    dfdlambda(1) = 0;
    dfdlambda(2) = 1;
    dfdlambda(3) = 1;
    break;
  case 1:
    dfdlambda(0) = 0;
    dfdlambda(1) = 0;
    dfdlambda(2) = -x[2];
    dfdlambda(3) = -x[3];
    break;
  case 2:
    dfdlambda(0) = -p[4] * pow(x[0], 3);
    dfdlambda(1) = -p[4] * pow(x[1], 3);
    dfdlambda(2) = 0;
    dfdlambda(3) = 0;
    break;
  case 3:
    dfdlambda(0) = ((1.0 / 2.0) * tanh(x[1]) + 1.0 / 2.0) * p[4];
    dfdlambda(1) = (-1.0 / 2.0 * tanh(x[0]) - 1.0 / 2.0) * p[4];
    dfdlambda(2) = 0;
    dfdlambda(3) = 0;
    break;
  case 4:
    dfdlambda(0) = ((1.0 / 2.0) * tanh(x[1]) + 1.0 / 2.0) * p[3] -
                   p[2] * pow(x[0], 3) + x[0] - x[2];
    dfdlambda(1) = -((1.0 / 2.0) * tanh(x[0]) + 1.0 / 2.0) * p[3] -
                   p[2] * pow(x[1], 3) + x[1] - x[3];
    dfdlambda(2) = 0;
    dfdlambda(3) = 0;
    break;
  }

  dfdxdx[0](0, 0) = -6 * p[2] * p[4] * x[0];
  dfdxdx[0](0, 1) = 0;
  dfdxdx[0](0, 2) = 0;
  dfdxdx[0](0, 3) = 0;
  dfdxdx[0](1, 0) =
      (1.0 / 2.0) * (2 - 2 * pow(tanh(x[0]), 2)) * tanh(x[0]) * p[3] * p[4];
  dfdxdx[0](1, 1) = 0;
  dfdxdx[0](1, 2) = 0;
  dfdxdx[0](1, 3) = 0;
  dfdxdx[0](2, 0) = 0;
  dfdxdx[0](2, 1) = 0;
  dfdxdx[0](2, 2) = 0;
  dfdxdx[0](2, 3) = 0;
  dfdxdx[0](3, 0) = 0;
  dfdxdx[0](3, 1) = 0;
  dfdxdx[0](3, 2) = 0;
  dfdxdx[0](3, 3) = 0;
  dfdxdx[1](0, 0) = 0;
  dfdxdx[1](0, 1) =
      -1.0 / 2.0 * (2 - 2 * pow(tanh(x[1]), 2)) * tanh(x[1]) * p[3] * p[4];
  dfdxdx[1](0, 2) = 0;
  dfdxdx[1](0, 3) = 0;
  dfdxdx[1](1, 0) = 0;
  dfdxdx[1](1, 1) = -6 * p[2] * p[4] * x[1];
  dfdxdx[1](1, 2) = 0;
  dfdxdx[1](1, 3) = 0;
  dfdxdx[1](2, 0) = 0;
  dfdxdx[1](2, 1) = 0;
  dfdxdx[1](2, 2) = 0;
  dfdxdx[1](2, 3) = 0;
  dfdxdx[1](3, 0) = 0;
  dfdxdx[1](3, 1) = 0;
  dfdxdx[1](3, 2) = 0;
  dfdxdx[1](3, 3) = 0;
  dfdxdx[2](0, 0) = 0;
  dfdxdx[2](0, 1) = 0;
  dfdxdx[2](0, 2) = 0;
  dfdxdx[2](0, 3) = 0;
  dfdxdx[2](1, 0) = 0;
  dfdxdx[2](1, 1) = 0;
  dfdxdx[2](1, 2) = 0;
  dfdxdx[2](1, 3) = 0;
  dfdxdx[2](2, 0) = 0;
  dfdxdx[2](2, 1) = 0;
  dfdxdx[2](2, 2) = 0;
  dfdxdx[2](2, 3) = 0;
  dfdxdx[2](3, 0) = 0;
  dfdxdx[2](3, 1) = 0;
  dfdxdx[2](3, 2) = 0;
  dfdxdx[2](3, 3) = 0;
  dfdxdx[3](0, 0) = 0;
  dfdxdx[3](0, 1) = 0;
  dfdxdx[3](0, 2) = 0;
  dfdxdx[3](0, 3) = 0;
  dfdxdx[3](1, 0) = 0;
  dfdxdx[3](1, 1) = 0;
  dfdxdx[3](1, 2) = 0;
  dfdxdx[3](1, 3) = 0;
  dfdxdx[3](2, 0) = 0;
  dfdxdx[3](2, 1) = 0;
  dfdxdx[3](2, 2) = 0;
  dfdxdx[3](2, 3) = 0;
  dfdxdx[3](3, 0) = 0;
  dfdxdx[3](3, 1) = 0;
  dfdxdx[3](3, 2) = 0;
  dfdxdx[3](3, 3) = 0;

  switch (var_param) {
  case 0:
    dfdxdlambda(0, 0) = 0;
    dfdxdlambda(0, 1) = 0;
    dfdxdlambda(0, 2) = 0;
    dfdxdlambda(0, 3) = 0;
    dfdxdlambda(1, 0) = 0;
    dfdxdlambda(1, 1) = 0;
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(1, 3) = 0;
    dfdxdlambda(2, 0) = 0;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = 0;
    dfdxdlambda(2, 3) = 0;
    dfdxdlambda(3, 0) = 0;
    dfdxdlambda(3, 1) = 0;
    dfdxdlambda(3, 2) = 0;
    dfdxdlambda(3, 3) = 0;
    break;
  case 1:
    dfdxdlambda(0, 0) = 0;
    dfdxdlambda(0, 1) = 0;
    dfdxdlambda(0, 2) = 0;
    dfdxdlambda(0, 3) = 0;
    dfdxdlambda(1, 0) = 0;
    dfdxdlambda(1, 1) = 0;
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(1, 3) = 0;
    dfdxdlambda(2, 0) = 0;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = -1;
    dfdxdlambda(2, 3) = 0;
    dfdxdlambda(3, 0) = 0;
    dfdxdlambda(3, 1) = 0;
    dfdxdlambda(3, 2) = 0;
    dfdxdlambda(3, 3) = -1;
    break;
  case 2:
    dfdxdlambda(0, 0) = -3 * p[4] * pow(x[0], 2);
    dfdxdlambda(0, 1) = 0;
    dfdxdlambda(0, 2) = 0;
    dfdxdlambda(0, 3) = 0;
    dfdxdlambda(1, 0) = 0;
    dfdxdlambda(1, 1) = -3 * p[4] * pow(x[1], 2);
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(1, 3) = 0;
    dfdxdlambda(2, 0) = 0;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = 0;
    dfdxdlambda(2, 3) = 0;
    dfdxdlambda(3, 0) = 0;
    dfdxdlambda(3, 1) = 0;
    dfdxdlambda(3, 2) = 0;
    dfdxdlambda(3, 3) = 0;
    break;
  case 3:
    dfdxdlambda(0, 0) = 0;
    dfdxdlambda(0, 1) = (1.0 / 2.0 - 1.0 / 2.0 * pow(tanh(x[1]), 2)) * p[4];
    dfdxdlambda(0, 2) = 0;
    dfdxdlambda(0, 3) = 0;
    dfdxdlambda(1, 0) = ((1.0 / 2.0) * pow(tanh(x[0]), 2) - 1.0 / 2.0) * p[4];
    dfdxdlambda(1, 1) = 0;
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(1, 3) = 0;
    dfdxdlambda(2, 0) = 0;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = 0;
    dfdxdlambda(2, 3) = 0;
    dfdxdlambda(3, 0) = 0;
    dfdxdlambda(3, 1) = 0;
    dfdxdlambda(3, 2) = 0;
    dfdxdlambda(3, 3) = 0;
    break;
  case 4:
    dfdxdlambda(0, 0) = -3 * p[2] * pow(x[0], 2) + 1;
    dfdxdlambda(0, 1) = (1.0 / 2.0 - 1.0 / 2.0 * pow(tanh(x[1]), 2)) * p[3];
    dfdxdlambda(0, 2) = -1;
    dfdxdlambda(0, 3) = 0;
    dfdxdlambda(1, 0) = ((1.0 / 2.0) * pow(tanh(x[0]), 2) - 1.0 / 2.0) * p[3];
    dfdxdlambda(1, 1) = -3 * p[2] * pow(x[1], 2) + 1;
    dfdxdlambda(1, 2) = 0;
    dfdxdlambda(1, 3) = -1;
    dfdxdlambda(2, 0) = 0;
    dfdxdlambda(2, 1) = 0;
    dfdxdlambda(2, 2) = 0;
    dfdxdlambda(2, 3) = 0;
    dfdxdlambda(3, 0) = 0;
    dfdxdlambda(3, 1) = 0;
    dfdxdlambda(3, 2) = 0;
    dfdxdlambda(3, 3) = 0;
    break;
  }

  if (mode != 4 && mode != 5 && mode != 6) { // no VE for eqp analysis
    unsigned int counter = xdim;
    // variational state (transform to matrix shape for easy producting)
    Eigen::MatrixXd state_dphidx = x(Eigen::seqN(counter, size_dphidx));
    state_dphidx.resize(xdim, xdim);
    counter += size_dphidx;
    Eigen::VectorXd state_dphidlambda =
        x(Eigen::seqN(counter, size_dphidlambda));
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

    // dphidlambda
    dxdt(Eigen::seqN(counter, xdim)) = dfdx * state_dphidlambda + dfdlambda;
    counter += size_dphidlambda;

    if (mode != 0 && !numerical_diff) {
      // dphidxdx
      for (int i = 0; i < xdim; i++) {
        dxdt(Eigen::seqN(counter, size_dphidx)) =
            dfdxdx[i] * state_dphidx.col(i).prod() * state_dphidx +
            dfdx * state_dphidxdx[i];
        counter += size_dphidx;
      }

      // dphidxdlambda
      for (int i = 0; i < xdim; i++) {
        dxdt(Eigen::seqN(counter, xdim)) =
            dfdxdx[i] * state_dphidlambda * state_dphidx +
            dfdxdlambda * state_dphidx.col(i) +
            dfdx * state_dphidxdlambda.col(i);
        counter += xdim;
      }
    }
  }
}