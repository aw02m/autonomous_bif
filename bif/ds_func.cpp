#include "ds_func.hpp"
#include "ds_derivatives.hpp"
#include "dynamical_system.hpp"
#include "eigensolver.hpp"
#include "runge_kutta.hpp"

void store_constant_state(dynamical_system &ds) {
  ds.dhdx = dhdx(ds);
  ds.dh_invdu = dh_invdu(ds);
  ds.dqdx = dqdx(ds);
}

void store_state(const Eigen::VectorXd &vp, dynamical_system &ds) {
  Eigen::VectorXd u = vp(Eigen::seqN(0, ds.udim));
  double tau = vp(ds.udim);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  Eigen::VectorXd O = Eigen::VectorXd::Zero(
      ds.size_dphidlambda + ds.size_dphidxdx + ds.size_dphidxdlambda);
  Eigen::VectorXd init_state(ds.xdim + ds.size_dphidx + ds.size_dphidlambda +
                             ds.size_dphidxdx + ds.size_dphidxdlambda);

  ds.u0 = u;
  ds.tauk(0) = tau;
  ds.xk[0] = h_inv(u, ds);

  I.resize(I.cols() * I.rows(), 1);
  init_state << ds.xk[0], I, O;
  Eigen::VectorXd sol;
  if (ds.use_classic_rk != false) {
    sol = integrate_rk45(variational_eq, 0, init_state, ds.tauk(0), ds);
  } else {
    sol = integrate(variational_eq, 0, init_state, ds.tauk(0), ds);
  }

  unsigned int counter = 0;
  ds.xk[1] = sol(Eigen::seqN(counter, ds.xdim));
  ds.fk[0] = f(0, ds.xk[1], ds);
  ds.dfdx[0] = dfdx(ds.xk[1], ds);
  counter += ds.xdim;

  Eigen::MatrixXd dphidx;
  dphidx = sol(Eigen::seqN(counter, ds.size_dphidx));
  dphidx.resize(ds.xdim, ds.xdim);
  ds.dphidx[0] = dphidx;
  counter += ds.size_dphidx;

  ds.dphidlambda[0] = sol(Eigen::seqN(counter, ds.size_dphidlambda));
  counter += ds.size_dphidlambda;

  Eigen::MatrixXd dphidxdxk;
  for (int i = 0; i < ds.xdim; i++) {
    dphidxdxk = sol(Eigen::seqN(counter + ds.size_dphidx * i, ds.size_dphidx));
    dphidxdxk.resize(ds.xdim, ds.xdim);
    ds.dphidxdx[0][i] = dphidxdxk;
    dphidxdxk.resize(ds.size_dphidx, 1);
  }
  counter += ds.size_dphidxdx;

  Eigen::MatrixXd dphidxdlambdak;
  dphidxdlambdak = sol(Eigen::seqN(counter, ds.size_dphidxdlambda));
  dphidxdlambdak.resize(ds.xdim, ds.xdim);
  ds.dphidxdlambda[0] = dphidxdlambdak;

  ds.dphidxdtau[0] = ds.dfdx[0] * ds.dphidx[0];

  ds.dTldu = dTldu(ds);
  ds.dTldtau = dTldtau(ds);
  ds.dTldlambda = dTldlambda(ds);
  ds.dqdu = dqdu(ds);
  ds.dqdtau = dqdtau(ds);
  ds.dqdlambda = dqdlambda(ds);
  ds.dTldudu = dTldudu(ds);
  ds.dTldudtau = dTldudtau(ds);
  ds.dTldudlambda = dTldudlambda(ds);

  ds.eigvals = eigenvalues(ds);

  ds.chara_poly = ds.dTldu - ds.mu * Eigen::MatrixXd::Identity(ds.udim, ds.udim);
}

Eigen::VectorXd variational_eq(double t, const Eigen::VectorXd &x,
                               const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim + ds.size_dphidx + ds.size_dphidlambda +
                      ds.size_dphidxdx + ds.size_dphidxdlambda);
  unsigned int counter = 0;

  Eigen::VectorXd state_x = x(Eigen::seqN(counter, ds.xdim));
  counter += ds.xdim;
  Eigen::MatrixXd state_dfdx = dfdx(state_x, ds);
  Eigen::VectorXd state_dfdlambda = dfdlambda(state_x, ds);
  std::vector<Eigen::MatrixXd> state_dfdxdx = dfdxdx(state_x, ds);
  Eigen::MatrixXd state_dfdxdlambda = dfdxdlambda(state_x, ds);

  Eigen::MatrixXd state_dphidx = x(Eigen::seqN(counter, ds.size_dphidx));
  state_dphidx.resize(ds.xdim, ds.xdim);
  counter += ds.size_dphidx;
  Eigen::VectorXd state_dphidlambda =
      x(Eigen::seqN(counter, ds.size_dphidlambda));
  counter += ds.size_dphidlambda;
  std::vector<Eigen::MatrixXd> state_dphidxdx(
      ds.xdim, Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));
  Eigen::MatrixXd temp;
  for (int i = 0; i < ds.xdim; i++) {
    temp = x(Eigen::seqN(counter + ds.size_dphidx * i, ds.size_dphidx));
    temp.resize(ds.xdim, ds.xdim);
    state_dphidxdx[i] = temp;
    temp.resize(ds.size_dphidx, 1);
  }
  counter += ds.size_dphidxdx;
  Eigen::MatrixXd state_dphidxdlambda =
      x(Eigen::seqN(counter, ds.size_dphidxdlambda));

  counter = 0;
  ret(Eigen::seqN(0, ds.xdim)) = f(t, state_x, ds);
  counter += ds.xdim;
  ret(Eigen::seqN(counter, ds.xdim * ds.xdim)) = state_dfdx * state_dphidx;
  counter += ds.size_dphidx;
  ret(Eigen::seqN(counter, ds.xdim)) =
      state_dfdx * state_dphidlambda + state_dfdlambda;
  counter += ds.size_dphidlambda;
  std::vector<Eigen::MatrixXd> DDF(ds.xdim,
                                   Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));
  for (int i = 0; i < ds.xdim; i++) {
    for (int j = 0; j < ds.xdim; j++) {
      DDF[i].col(j) = state_dfdxdx[j] * state_dphidx.col(i);
    }
    // std::cout << state_dfdx * state_dphidxdx[i] << std::endl;
    ret(Eigen::seqN(counter + ds.size_dphidx * i, ds.size_dphidx)) =
        state_dfdx * state_dphidxdx[i] + DDF[i] * state_dphidx;
    // std::cout << state_dfdx * state_dphidxdx[i] + DDF[i] * state_dphidx
    //           << std::endl;
    // std::cout << ret(Eigen::seqN(counter + ds.size_dphidx * i,
    // ds.size_dphidx)); exit(0);
  }
  counter += ds.size_dphidxdx;
  for (int i = 0; i < ds.xdim; i++) {
    ret(Eigen::seqN(counter + ds.xdim * i, ds.xdim)) =
        state_dfdx * state_dphidxdlambda.col(i) + DDF[i] * state_dphidlambda +
        state_dfdxdlambda.col(i);
  }

  return ret;
}

Eigen::VectorXd func_newton(const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.udim + ds.period + 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);

  ret(Eigen::seqN(0, ds.udim)) = h(ds.xk[ds.period], ds) - ds.u0;
  ret(ds.udim) = q(ds.xk[ds.period], ds);
  ret(ds.udim + 1) = ds.chara_poly.determinant();

  return ret;
}

double det_derivative(const Eigen::MatrixXd &A, const Eigen::MatrixXd &dA,
                      const dynamical_system &ds) {
  Eigen::MatrixXd temp(ds.udim, ds.udim);
  double ret = 0;

  for (int i = 0; i < ds.udim; i++) {
    temp = A;
    temp.col(i) = dA.col(i);
    ret += temp.determinant();
  }

  return ret;
}

Eigen::MatrixXd jac_newton(const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.udim + ds.period + 1, ds.udim + ds.period + 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.udim, ds.udim);

  ret(Eigen::seqN(0, ds.udim), Eigen::seqN(0, ds.udim)) = ds.dTldu - I;
  ret(Eigen::seqN(0, ds.udim), ds.udim) = ds.dTldtau;
  ret(Eigen::seqN(0, ds.udim), ds.udim + 1) = ds.dTldlambda;
  ret(ds.udim, Eigen::seqN(0, ds.udim)) = ds.dqdu;
  ret(ds.udim, ds.udim) = ds.dqdtau(0, 0);
  ret(ds.udim, ds.udim + 1) = ds.dqdlambda(0, 0);
  ret(ds.udim + 1, Eigen::seqN(0, ds.udim)) = dchidu(ds);
  ret(ds.udim + 1, ds.udim) = dchidtau(ds);
  ret(ds.udim + 1, ds.udim + 1) = dchidlambda(ds);

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

// Eigen::VectorXd variational_init1(double t, const Eigen::MatrixXd &dphidx,
//                                   const Eigen::MatrixXd &dfdx,
//                                   const dynamical_system &ds) {
//   Eigen::MatrixXd ret = dfdx * dphidx;
//   ret.resize(ds.xdim * ds.xdim, 1);
//   return ret;
// }

// Eigen::VectorXd variational_para1(double t, const Eigen::VectorXd
// &dphidlambda,
//                                   const Eigen::VectorXd &x,
//                                   const Eigen::MatrixXd &dfdx,
//                                   const dynamical_system &ds) {
//   // Eigen::VectorXd ret;
//   // ret = dfdx(phi, ds) * x + dfdlambda(phi, ds);
//   // return ret;
//   return dfdx * x + dfdlambda(x, ds);
// }

Eigen::MatrixXd dTldu(const dynamical_system &ds) {
  return ds.dhdx * ds.dphidx[0] * ds.dh_invdu;
}

Eigen::VectorXd dTldtau(const dynamical_system &ds) {
  return ds.dhdx * ds.fk[0];
}

Eigen::VectorXd dTldlambda(const dynamical_system &ds) {
  return ds.dhdx * ds.dphidlambda[0];
}

Eigen::MatrixXd dqdu(const dynamical_system &ds) {
  return ds.dqdx * ds.dphidx[0] * ds.dh_invdu;
}

Eigen::MatrixXd dqdtau(const dynamical_system &ds) {
  return ds.dqdx * ds.fk[0];
}

Eigen::VectorXd dqdlambda(const dynamical_system &ds) {
  return ds.dqdx * ds.dphidlambda[0];
}

std::vector<Eigen::MatrixXd> dTldudu(const dynamical_system &ds) {
  std::vector<Eigen::MatrixXd> ret = std::vector<Eigen::MatrixXd>(
      ds.udim, Eigen::MatrixXd::Zero(ds.udim, ds.udim));
  for (int i = 0; i < ds.udim; i++) {
    if (i != ds.p_index) {
      if (i < ds.p_index) {
        ret[i] = ds.dhdx * ds.dphidxdx[0][i] * ds.dh_invdu;
      } else {
        ret[i - 1] = ds.dhdx * ds.dphidxdx[0][i] * ds.dh_invdu;
      }
    }
  }
  return ret;
}

Eigen::MatrixXd dTldudtau(const dynamical_system &ds) {
  return ds.dhdx * ds.dphidxdtau[0] * ds.dh_invdu;
}

Eigen::MatrixXd dTldudlambda(const dynamical_system &ds) {
  return ds.dhdx * ds.dphidxdlambda[0] * ds.dh_invdu;
}

Eigen::MatrixXd dchidu(const dynamical_system &ds) {
  Eigen::MatrixXd ret(1, ds.udim);
  for (int i = 0; i < ds.udim; i++) {
    ret(0, i) = det_derivative(ds.chara_poly, ds.dTldudu[i], ds);
  }
  return ret;
}

double dchidtau(const dynamical_system &ds) {
  return det_derivative(ds.chara_poly, ds.dTldudtau, ds);
}

double dchidlambda(const dynamical_system &ds) {
  return det_derivative(ds.chara_poly, ds.dTldudlambda, ds);
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