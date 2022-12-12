#include "dynamical_system.hpp"
#include "essential.hpp"

std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
dynamical_system::newton_FJ(const Eigen::VectorXd &v) {
  store_constant_state();
  switch (mode) {
  case 3:
    if (numerical_diff) {
      store_states_numeric(v);
    } else {
      store_states(v);
    }
    return std::make_tuple(newton_F_NS(), newton_J_NS());
  case 1:
    if (numerical_diff) {
      store_states_numeric(v);
    } else {
      store_states(v);
    }
    return std::make_tuple(newton_F_G(), newton_J_G());
  case 2:
    if (numerical_diff) {
      store_states_numeric(v);
    } else {
      store_states(v);
    }
    return std::make_tuple(newton_F_PD(), newton_J_PD());
  case 0:
    store_states_fix(v);
    return std::make_tuple(newton_F_fix(), newton_J_fix());
  case 6:
    store_states_eqp(v);
    return std::make_tuple(newton_F_eqp_H(), newton_J_eqp_H());
  case 5:
    store_states_eqp(v);
    return std::make_tuple(newton_F_eqp_G(), newton_J_eqp_G());
  case 4:
    store_states_eqp(v);
    return std::make_tuple(newton_F_eqp(), newton_J_eqp());
  }
  // exception
  return std::make_tuple(Eigen::VectorXd::Zero(xdim),
                         Eigen::MatrixXd::Zero(xdim, xdim));
}

Eigen::VectorXd dynamical_system::newton_F_G() {
  Eigen::VectorXd ret(xdim + 2);

  ret(Eigen::seqN(0, xdim)) = xk[1] - xk[0];
  ret(xdim) = q(xk[1]);
  ret(xdim + 1) =
      det_derivative(chara_poly, -Eigen::MatrixXd::Identity(xdim, xdim), xdim);

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_G() {
  Eigen::MatrixXd ret(xdim + 2, xdim + 2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTdx - I;
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(Eigen::seqN(0, xdim), xdim + 1) = dTdlambda;
  ret(xdim, Eigen::seqN(0, xdim)) = dpidx;
  ret(xdim, xdim) = dpidtau(0, 0);
  ret(xdim, xdim + 1) = dpidlambda(0, 0);

  Eigen::MatrixXd temp(xdim, xdim);
  Eigen::MatrixXd dtemp(xdim, xdim);
  // compute dchi/dmudx caz bif cond. is dchi/dmu
  for (int i = 0; i < xdim; i++) {
    ret(xdim + 1, i) = 0;
    for (int j = 0; j < xdim; j++) {
      temp = chara_poly;
      temp.col(j) = Eigen::VectorXd::Zero(xdim);
      temp(j, j) = -1;
      dtemp = dTdxdx[i];
      dtemp.col(j) = Eigen::VectorXd::Zero(xdim);
      ret(xdim + 1, i) += det_derivative(temp, dtemp, xdim);
    }
  }
  // dchi/dmudtau
  ret(xdim + 1, xdim) = 0;
  for (int j = 0; j < xdim; j++) {
    temp = chara_poly;
    temp.col(j) = Eigen::VectorXd::Zero(xdim);
    temp(j, j) = -1;
    dtemp = dTdxdtau;
    dtemp.col(j) = Eigen::VectorXd::Zero(xdim);
    ret(xdim + 1, xdim) += det_derivative(temp, dtemp, xdim);
  }
  // dchi/dmudlambda
  ret(xdim + 1, xdim + 1) = 0;
  for (int j = 0; j < xdim; j++) {
    temp = chara_poly;
    temp.col(j) = Eigen::VectorXd::Zero(xdim);
    temp(j, j) = -1;
    dtemp = dTdxdlambda;
    dtemp.col(j) = Eigen::VectorXd::Zero(xdim);
    ret(xdim + 1, xdim + 1) += det_derivative(temp, dtemp, xdim);
  }

  return ret;
}

Eigen::VectorXd dynamical_system::newton_F_PD() {
  Eigen::VectorXd ret(xdim + 2);

  ret(Eigen::seqN(0, xdim)) = xk[1] - xk[0];
  ret(xdim) = q(xk[1]);
  ret(xdim + 1) = chara_poly.determinant();

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_PD() {
  Eigen::MatrixXd ret(xdim + 2, xdim + 2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTdx - I;
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(Eigen::seqN(0, xdim), xdim + 1) = dTdlambda;
  ret(xdim, Eigen::seqN(0, xdim)) = dpidx;
  ret(xdim, xdim) = dpidtau(0, 0);
  ret(xdim, xdim + 1) = dpidlambda(0, 0);

  Eigen::MatrixXd dchidx = Eigen::MatrixXd::Zero(1, xdim);
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dTdxdx[i], xdim);
  }
  double dchidtau = det_derivative(chara_poly, dTdxdtau, xdim);
  double dchidlambda = det_derivative(chara_poly, dTdxdlambda, xdim);
  ret(xdim + 1, Eigen::seqN(0, xdim)) = dchidx;
  ret(xdim + 1, xdim) = dchidtau;
  ret(xdim + 1, xdim + 1) = dchidlambda;

  return ret;
}

Eigen::VectorXd dynamical_system::newton_F_NS() {
  Eigen::VectorXd ret(xdim + 2);

  ret(Eigen::seqN(0, xdim)) = xk[period] - xk[0];
  ret(xdim) = q(xk[period]);
  ret(xdim + 1) = chara_poly.determinant();

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_NS() {
  Eigen::MatrixXd ret(xdim + 2, xdim + 2);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) =
      dTdx - Eigen::MatrixXd::Identity(xdim, xdim);
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(Eigen::seqN(0, xdim), xdim + 1) = dTdlambda;
  ret(xdim, Eigen::seqN(0, xdim)) = dpidx;
  ret(xdim, xdim) = dpidtau(0, 0);
  ret(xdim, xdim + 1) = dpidlambda(0, 0);

  Eigen::MatrixXd dchidx = Eigen::MatrixXd::Zero(1, xdim);
  Eigen::MatrixXd temp(bialt_dim, bialt_dim);
  for (int i = 0; i < xdim; i++) {
    temp = bialt_prod_square_derivative(dTdx, dTdxdx[i], xdim, bialt_dim);
    dchidx(0, i) = det_derivative(chara_poly, temp, bialt_dim);
  }
  temp = bialt_prod_square_derivative(dTdx, dTdxdtau, xdim, bialt_dim);
  double dchidtau = det_derivative(chara_poly, temp, bialt_dim);
  temp = bialt_prod_square_derivative(dTdx, dTdxdlambda, xdim, bialt_dim);
  double dchidlambda = det_derivative(chara_poly, temp, bialt_dim);
  ret(xdim + 1, Eigen::seqN(0, xdim)) = dchidx;
  ret(xdim + 1, xdim) = dchidtau;
  ret(xdim + 1, xdim + 1) = dchidlambda;

  return ret;
}

Eigen::VectorXd dynamical_system::newton_F_fix() {
  Eigen::VectorXd ret(xdim + 1);

  ret(Eigen::seqN(0, xdim)) = xk[1] - xk[0];
  ret(xdim) = q(xk[1]);

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_fix() {
  Eigen::MatrixXd ret(xdim + 1, xdim + 1);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) =
      dTdx - Eigen::MatrixXd::Identity(xdim, xdim);
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(xdim, Eigen::seqN(0, xdim)) = dpidx;
  ret(xdim, xdim) = dpidtau(0, 0);

  return ret;
}

Eigen::VectorXd dynamical_system::newton_F_eqp() {
  Eigen::VectorXd ret(xdim);

  ret(Eigen::seqN(0, xdim)) = f;

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_eqp() {
  Eigen::MatrixXd ret(xdim, xdim);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dfdx;

  return ret;
}

Eigen::VectorXd dynamical_system::newton_F_eqp_G() {
  Eigen::VectorXd ret(xdim + 1);

  ret(Eigen::seqN(0, xdim)) = f;
  ret(xdim) = chara_poly.determinant();

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_eqp_G() {
  Eigen::MatrixXd ret(xdim + 1, xdim + 1);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dfdx;
  ret(Eigen::seqN(0, xdim), xdim) = dfdlambda;

  Eigen::MatrixXd dchidx = Eigen::MatrixXd::Zero(1, xdim);
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dfdxdx[i], bialt_dim);
  }
  double dchidlambda = det_derivative(chara_poly, dfdxdlambda, bialt_dim);
  ret(xdim, Eigen::seqN(0, xdim)) = dchidx;
  ret(xdim, xdim) = dchidlambda;

  return ret;
}

Eigen::VectorXd dynamical_system::newton_F_eqp_H() {
  Eigen::VectorXd ret(xdim + 1);

  ret(Eigen::seqN(0, xdim)) = f;
  ret(xdim) = chara_poly.determinant();

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_eqp_H() {
  Eigen::MatrixXd ret(xdim + 1, xdim + 1);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dfdx;
  ret(Eigen::seqN(0, xdim), xdim) = dfdlambda;

  Eigen::MatrixXd dchidx = Eigen::MatrixXd::Zero(1, xdim);
  Eigen::MatrixXd temp(bialt_dim, bialt_dim);
  for (int i = 0; i < xdim; i++) {
    temp = biproduct_derivative(dfdxdx[i], xdim, bialt_dim);
    dchidx(0, i) = det_derivative(chara_poly, temp, bialt_dim);
  }
  temp = biproduct_derivative(dfdxdlambda, xdim, bialt_dim);
  double dchidlambda = det_derivative(chara_poly, temp, bialt_dim);
  ret(xdim, Eigen::seqN(0, xdim)) = dchidx;
  ret(xdim, xdim) = dchidlambda;

  return ret;
}