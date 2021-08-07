#include "ds_func.hpp"

void store_constant_state(dynamical_system &ds) {
  ds.dqdx = dqdx(ds);
  ds.dhdx = dhdx(ds);
  ds.dh_invdu = dh_invdu(ds);
}

void store_state(const Eigen::VectorXd &vp, dynamical_system &ds) {
  Eigen::VectorXd u = vp(Eigen::seqN(0, ds.udim));
  double tau = vp(ds.udim);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  Eigen::VectorXd init_state(ds.xdim + ds.xdim * ds.xdim);

  ds.u0 = u;
  ds.tauk(0) = tau;
  ds.xk[0] = h_inv(u, ds);

  I.resize(I.cols() * I.rows(), 1);
  init_state << ds.xk[0], I;
  Eigen::VectorXd sol;
  if (ds.use_classic_rk != false) {
    sol = integrate_rk45(variational_eq, 0, init_state, ds.tauk(0), ds);
  } else {
    sol = integrate(variational_eq, 0, init_state, ds.tauk(0), ds);
  }

  ds.xk[1] = sol(Eigen::seqN(0, ds.xdim));

  Eigen::MatrixXd dphidx;
  dphidx = sol(Eigen::seqN(ds.xdim, ds.xdim * ds.xdim));
  dphidx.resize(ds.xdim, ds.xdim);
  ds.dphidx[0] = dphidx;

  ds.fk[0] = f(0, ds.xk[1], ds);
  ds.dfdx[0] = dfdx(ds.xk[1], ds);

  ds.dTldu = dTldu(ds);
}

Eigen::VectorXd func_newton(const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.udim + ds.period);

  ret(Eigen::seqN(0, ds.udim)) = ds.u0 - h(ds.xk[ds.period], ds);
  ret(ds.udim) = q(ds.xk[ds.period], ds);

  return ret;
}

Eigen::MatrixXd jac_newton(const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.udim + ds.period, ds.udim + ds.period);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.udim, ds.udim);

  ret(Eigen::seqN(0, ds.udim), Eigen::seqN(0, ds.udim)) = I - ds.dTldu;
  ret(Eigen::seqN(0, ds.udim), ds.udim) = dTldtau(ds);
  ret(ds.udim, Eigen::seqN(0, ds.udim)) = dqdu(ds);
  ret(ds.udim, ds.udim) = dqdtau(ds)(0, 0);

  return ret;
}

Eigen::VectorXd T(const Eigen::VectorXd &x, double tau,
                  const dynamical_system &ds) {
  return integrate(f, 0, x, tau, ds);
}

double q(const Eigen::VectorXd &x, const dynamical_system &ds) {
  return x(ds.p_index) - ds.p_place;
}

Eigen::VectorXd h(const Eigen::VectorXd &x, const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.udim);

  for (int i = 0; i < ds.xdim; i++) {
    if (i != ds.p_index) {
      if (i < ds.p_index) {
        ret(i) = x(i);
      } else {
        ret(i - 1) = x(i);
      }
    }
  }

  return ret;
}

Eigen::VectorXd h_inv(const Eigen::VectorXd &u, const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim);

  for (int i = 0; i < ds.xdim; i++) {
    if (i != ds.p_index) {
      if (i < ds.p_index) {
        ret(i) = u(i);
      } else {
        ret(i) = u(i - 1);
      }
    } else {
      ret(i) = ds.p_place;
    }
  }

  return ret;
}

Eigen::MatrixXd dqdx(const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(1, ds.xdim);
  ret(0, ds.p_index) = 1.0;
  return ret;
}

Eigen::MatrixXd dhdx(const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  removeRow(ret, ds.p_index);
  return ret;
}

Eigen::MatrixXd dh_invdu(const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  removeCol(ret, ds.p_index);
  return ret;
}

Eigen::VectorXd variational_eq(double t, const Eigen::VectorXd &x,
                               const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim + ds.xdim * ds.xdim);

  Eigen::VectorXd state_x = x(Eigen::seqN(0, ds.xdim));

  Eigen::MatrixXd state_dphidx(ds.xdim, ds.xdim);
  state_dphidx = x(Eigen::seqN(ds.xdim, ds.xdim * ds.xdim));
  state_dphidx.resize(ds.xdim, ds.xdim);

  ret(Eigen::seqN(0, ds.xdim)) = f(t, state_x, ds);
  ret(Eigen::seqN(ds.xdim, ds.xdim * ds.xdim)) =
      variational_init1(t, state_x, state_dphidx, ds);

  return ret;
}

Eigen::VectorXd variational_init1(double t, const Eigen::VectorXd &phi,
                                  const Eigen::MatrixXd &dphidx,
                                  const dynamical_system &ds) {
  Eigen::MatrixXd ret = dfdx(phi, ds) * dphidx;
  ret.resize(ds.xdim * ds.xdim, 1);
  return ret;
}

Eigen::MatrixXd dTldu(const dynamical_system &ds) {
  return ds.dhdx * ds.dphidx[0] * ds.dh_invdu;
}

Eigen::VectorXd dTldtau(const dynamical_system &ds) {
  return ds.dhdx * ds.fk[0];
}

Eigen::MatrixXd dqdu(const dynamical_system &ds) {
  return ds.dqdx * ds.dphidx[0] * ds.dh_invdu;
}

Eigen::MatrixXd dqdtau(const dynamical_system &ds) {
  return ds.dqdx * ds.fk[0];
}

void removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove) {
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

  matrix.conservativeResize(numRows, numCols);
}

void removeCol(Eigen::MatrixXd &matrix, unsigned int colToRemove) {
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
        matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}

//////////////////////////////////////////////////////////////////////
///// Implementation example of projection h & h_inv using Eigen /////
//////////////////////////////////////////////////////////////////////

// // Easy
// Eigen::VectorXd h(const Eigen::VectorXd& x, const dynamical_system& ds)
// {
//   Eigen::VectorXd ret(ds.udim);
//   ret << x(0), x(2);
//   return ret;
// }

// Eigen::VectorXd h_inv(const Eigen::VectorXd& u, const dynamical_system& ds)
// {
//   Eigen::VectorXd ret(ds.xdim);
//   ret << u(0), ds.p_place, u(1);
//   return ret;
// }

//// Rich but slower
// Eigen::VectorXd h(Eigen::VectorXd &x, dynamical_system &ds) {
//   unsigned int num_rows = x.rows() - 1;
//   Eigen::VectorXd ret(x.rows() - 1);

//   if (ds.p_index < num_rows) {
//     ret.segment(ds.p_index, num_rows - ds.p_index) =
//         ret.bottomRows(num_rows - ds.p_index);
//   } else {
//     std::cerr << "psec_index is out of range in the function h." <<
//     std::endl; std::exit(1);
//   }

//   ret.conservativeResize(num_rows);
//   return ret;
// }

// Eigen::VectorXd h_inv(Eigen::VectorXd &u, dynamical_system &ds) {
//   unsigned int num_rows = u.rows() + 1;
//   Eigen::VectorXd ret(u.rows() + 1);
//   ret.head(u.rows()) = u;

//   if (ds.p_index < num_rows) {
//     ret.segment(ds.p_index + 1, num_rows - ds.p_index - 1) =
//         u.bottomRows(u.rows() - ds.p_index);
//     ret(ds.p_index) = ds.p_place;
//   } else {
//     std::cerr << "psec_index is out of range in the function h_inv."
//               << std::endl;
//     std::exit(1);
//   }

//   ret.conservativeResize(num_rows);
//   return ret;
// }