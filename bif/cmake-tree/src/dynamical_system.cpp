#include "dynamical_system.hpp"
#include "essential.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["x0"].size();
  bialt_dim = 0;
  for (int i = 0; i < xdim; i++)
    bialt_dim += i;

  p_index = json["p_index"];
  p_place = json["p_place"];

  use_classic_rk = json["use_classic_rk"];
  rk_h = json["rk_h"];
  rk_atol = json["rk_atol"];
  rk_rtol = json["rk_rtol"];

  period = json["period"];
  inc_param = json["inc_param"];
  var_param = json["var_param"];
  delta_inc = json["delta_inc"];
  inc_iter = json["inc_iter"];
  max_iter = json["max_iter"];
  numerical_diff = json["numerical_diff"];
  diff_strip = json["diff_strip"];
  eps = json["eps"];
  explode = json["explode"];

  mode = json["mode"];

  out_path = json["out_path"];
  json_out_path = json["json_out_path"];

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["x0"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;

  tau = json["tau"];

  std::vector<double> params_arr = json["params"];
  Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
  this->p = params;

  xk = std::vector<Eigen::VectorXd>(period + 1, Eigen::VectorXd::Zero(xdim));
  f = Eigen::VectorXd::Zero(xdim);
  dfdx = Eigen::MatrixXd::Zero(xdim, xdim);
  dfdlambda = Eigen::VectorXd::Zero(xdim);
  dfdxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));
  dfdxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);

  dphidx = Eigen::MatrixXd::Zero(xdim, xdim);
  dphidlambda = Eigen::VectorXd::Zero(xdim);
  dphidxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));
  dphidxdtau = Eigen::MatrixXd::Zero(xdim, xdim);
  dphidxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);

  size_dphidx = xdim * xdim;
  size_dphidlambda = xdim;
  size_dphidxdx = xdim * xdim * xdim;
  size_dphidxdlambda = xdim * xdim;
}

double dynamical_system::q(const Eigen::VectorXd &x) {
  return x(p_index) - p_place;
}

void dynamical_system::store_states(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd O = Eigen::VectorXd::Zero(size_dphidlambda + size_dphidxdx +
                                            size_dphidxdlambda);
  Eigen::VectorXd state(xdim + size_dphidx + size_dphidlambda + size_dphidxdx +
                        size_dphidxdlambda);

  xk[0] = v(Eigen::seqN(0, xdim));
  tau = v(xdim);
  p(var_param) = v(xdim + 1);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I, O;

  if (!use_classic_rk) {
    exit(1);
  } else {
    rk45_classic(state, 0, tau);
  }

  unsigned int counter = 0;
  xk[1] = state(Eigen::seqN(counter, xdim));
  counter += xdim;

  Eigen::MatrixXd dphidx_buf;
  dphidx_buf = state(Eigen::seqN(counter, size_dphidx));
  dphidx_buf.resize(xdim, xdim);
  dphidx = dphidx_buf;
  counter += size_dphidx;

  dphidlambda = state(Eigen::seqN(counter, size_dphidlambda));
  counter += size_dphidlambda;

  Eigen::MatrixXd dphidxdx_buf;
  for (int i = 0; i < xdim; i++) {
    dphidxdx_buf = state(Eigen::seqN(counter + size_dphidx * i, size_dphidx));
    dphidxdx_buf.resize(xdim, xdim);
    dphidxdx[i] = dphidxdx_buf;
    dphidxdx_buf.resize(size_dphidx, 1);
  }
  counter += size_dphidxdx;

  Eigen::MatrixXd dphidxdlambda_buf;
  dphidxdlambda_buf = state(Eigen::seqN(counter, size_dphidxdlambda));
  dphidxdlambda_buf.resize(xdim, xdim);
  dphidxdlambda = dphidxdlambda_buf;
  counter += size_dphidxdlambda;

  dphidxdtau = dfdx * dphidx;

  dTdx = dphidx;
  dTdtau = f;
  dTdlambda = dphidlambda;
  dpidx = dqdx * dphidx;
  dpidtau = dqdx * f;
  dpidlambda = dqdx * dphidlambda;
  dTdxdx = dphidxdx;
  dTdxdtau = dphidxdtau;
  dTdxdlambda = dphidxdlambda;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTdx).eigenvalues();
  theta = std::arg(eigvals(nearest_idx(eigvals, 0)));

  switch (mode) {
  case 1:
    chara_poly = dTdx - Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 2:
    chara_poly = dTdx + Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 3:
    chara_poly = bialt_prod_square(dTdx, xdim, bialt_dim) -
                 Eigen::MatrixXd::Identity(bialt_dim, bialt_dim);
    break;
  }
}

void dynamical_system::store_states_numeric(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd O = Eigen::VectorXd::Zero(size_dphidlambda + size_dphidxdx +
                                            size_dphidxdlambda);
  Eigen::VectorXd state(xdim + size_dphidx + size_dphidlambda + size_dphidxdx +
                        size_dphidxdlambda);
  Eigen::VectorXd state_init(xdim + size_dphidx + size_dphidlambda +
                             size_dphidxdx + size_dphidxdlambda);
  Eigen::VectorXd h = Eigen::VectorXd::Zero(xdim);
  Eigen::MatrixXd Ah;

  xk[0] = v(Eigen::seqN(0, xdim));
  tau = v(xdim);
  p(var_param) = v(xdim + 1);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I, O;
  state_init = state;

  if (!use_classic_rk) {
    exit(1);
  } else {
    rk45_classic(state, 0, tau);
  }

  unsigned int counter = 0;
  xk[1] = state(Eigen::seqN(counter, xdim));
  counter += xdim;

  Eigen::MatrixXd dphidx_buf;
  dphidx_buf = state(Eigen::seqN(counter, size_dphidx));
  dphidx_buf.resize(xdim, xdim);
  dphidx = dphidx_buf;
  counter += size_dphidx;

  dphidlambda = state(Eigen::seqN(counter, size_dphidlambda));
  counter += size_dphidlambda;

  // dTdxdx
  for (int i = 0; i < xdim; i++) {
    state = state_init;
    state(i) += diff_strip;
    if (!use_classic_rk) {
      exit(1);
    } else {
      rk45_classic(state, 0, tau);
    }
    Ah = state(Eigen::seqN(xdim, size_dphidx));
    Ah.resize(xdim, xdim);
    dphidxdx[i] = (Ah - dphidx) / diff_strip;
    Ah.resize(size_dphidx, 1);
  }

  // dTdxdlambda
  static double lambda_temp;
  lambda_temp = p(var_param);
  p(var_param) += diff_strip;
  state = state_init;
  if (!use_classic_rk) {
    exit(1);
  } else {
    rk45_classic(state, 0, tau);
  }
  Ah = state(Eigen::seqN(xdim, size_dphidx));
  Ah.resize(xdim, xdim);
  dphidxdlambda = (Ah - dphidx) / diff_strip;
  p(var_param) = lambda_temp;

  dphidxdtau = dfdx * dphidx;

  dTdx = dphidx;
  dTdtau = f;
  dTdlambda = dphidlambda;
  dpidx = dqdx * dphidx;
  dpidtau = dqdx * f;
  dpidlambda = dqdx * dphidlambda;
  dTdxdx = dphidxdx;
  dTdxdtau = dphidxdtau;
  dTdxdlambda = dphidxdlambda;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTdx).eigenvalues();
  theta = std::arg(eigvals(nearest_idx(eigvals, 0)));

  switch (mode) {
  case 1:
    chara_poly = dTdx - Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 2:
    chara_poly = dTdx + Eigen::MatrixXd::Identity(xdim, xdim);
    break;
  case 3:
    chara_poly = bialt_prod_square(dTdx, xdim, bialt_dim) -
                 Eigen::MatrixXd::Identity(bialt_dim, bialt_dim);
    break;
  }
}

void dynamical_system::store_constant_state() {
  dqdx = Eigen::MatrixXd::Zero(1, xdim);
  dqdx(0, p_index) = 1.0;
}

Eigen::VectorXd dynamical_system::newton_F() {
  Eigen::VectorXd ret(xdim + 2);

  ret(Eigen::seqN(0, xdim)) = xk[1] - xk[0];
  ret(xdim) = q(xk[1]);
  ret(xdim + 1) = chara_poly.determinant();

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J() {
  Eigen::MatrixXd ret(xdim + 2, xdim + 2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::MatrixXd dchidx = Eigen::MatrixXd::Zero(1, xdim);
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dTdxdx[i], xdim);
  }
  double dchidtau = det_derivative(chara_poly, dTdxdtau, xdim);
  double dchidlambda = det_derivative(chara_poly, dTdxdlambda, xdim);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTdx - I;
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(Eigen::seqN(0, xdim), xdim + 1) = dTdlambda;
  ret(xdim, Eigen::seqN(0, xdim)) = dpidx;
  ret(xdim, xdim) = dpidtau(0, 0);
  ret(xdim, xdim + 1) = dpidlambda(0, 0);
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

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) =
      dTdx - Eigen::MatrixXd::Identity(xdim, xdim);
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(Eigen::seqN(0, xdim), xdim + 1) = dTdlambda;
  ret(xdim, Eigen::seqN(0, xdim)) = dpidx;
  ret(xdim, xdim) = dpidtau(0, 0);
  ret(xdim, xdim + 1) = dpidlambda(0, 0);
  ret(xdim + 1, Eigen::seqN(0, xdim)) = dchidx;
  ret(xdim + 1, xdim) = dchidtau;
  ret(xdim + 1, xdim + 1) = dchidlambda;

  return ret;
}

void dynamical_system::store_states_fix(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd state(xdim + size_dphidx);

  xk[0] = v(Eigen::seqN(0, xdim));
  tau = v(xdim);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I;

  if (!use_classic_rk) {
    exit(1);
  } else {
    rk45_classic(state, 0, tau);
  }

  unsigned int counter = 0;
  xk[1] = state(Eigen::seqN(counter, xdim));
  counter += xdim;

  Eigen::MatrixXd dphidx_buf;
  dphidx_buf = state(Eigen::seqN(counter, size_dphidx));
  dphidx_buf.resize(xdim, xdim);
  dphidx = dphidx_buf;
  counter += size_dphidx;

  dTdx = dphidx;
  dTdtau = f;
  dpidx = dqdx * dphidx;
  dpidtau = dqdx * f;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTdx).eigenvalues();
  theta = std::arg(eigvals(nearest_idx(eigvals, 0)));
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

void dynamical_system::store_states_eqp(const Eigen::VectorXd &v) {
  xk[0] = v(Eigen::seqN(0, xdim));
  if (mode != 4) {
    p(var_param) = v(xdim);
  }

  Eigen::VectorXd dummy = Eigen::VectorXd::Zero(xdim);
  double t = 0;
  sys_func(xk[0], t);

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dfdx).eigenvalues();
  theta = std::arg(eigvals(nearest_idx(eigvals, 0)));

  switch (mode) {
  case 4:
    break;
  case 5:
    chara_poly = dfdx;
    break;
  case 6:
    chara_poly = biproduct(dfdx, xdim, bialt_dim);
    break;
  }
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

Eigen::VectorXd dynamical_system::newton_F_eqp_H() {
  Eigen::VectorXd ret(xdim + 1);

  ret(Eigen::seqN(0, xdim)) = f;
  ret(xdim) = chara_poly.determinant();

  // debug(ret);

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_eqp_H() {
  Eigen::MatrixXd ret(xdim + 1, xdim + 1);
  Eigen::MatrixXd dchidx = Eigen::MatrixXd::Zero(1, xdim);
  Eigen::MatrixXd temp(bialt_dim, bialt_dim);
  for (int i = 0; i < xdim; i++) {
    temp = biproduct_derivative(dfdxdx[i], xdim, bialt_dim);
    dchidx(0, i) = det_derivative(chara_poly, temp, bialt_dim);
  }
  temp = biproduct_derivative(dfdxdlambda, xdim, bialt_dim);
  double dchidlambda = det_derivative(chara_poly, temp, bialt_dim);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dfdx;
  ret(Eigen::seqN(0, xdim), xdim) = dfdlambda;
  ret(xdim, Eigen::seqN(0, xdim)) = dchidx;
  ret(xdim, xdim) = dchidlambda;

  return ret;
}

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
  case 2:
    if (numerical_diff) {
      store_states_numeric(v);
    } else {
      store_states(v);
    }
    return std::make_tuple(newton_F(), newton_J());
  case 0:
    store_states_fix(v);
    return std::make_tuple(newton_F_fix(), newton_J_fix());
  case 6:
    store_states_eqp(v);
    return std::make_tuple(newton_F_eqp_H(), newton_J_eqp_H());
  // case 5:
  //   store_states_eqp(v);
  //   return std::make_tuple(newton_F_eqp_G(), newton_J_eqp_G());
  case 4:
    store_states_eqp(v);
    return std::make_tuple(newton_F_eqp(), newton_J_eqp());
  }
  // exception
  return std::make_tuple(Eigen::VectorXd::Zero(xdim),
                         Eigen::MatrixXd::Zero(xdim, xdim));
}

void dynamical_system::rk45_classic(Eigen::VectorXd &x, double t0,
                                    double t_end) {
  Eigen::VectorXd k1 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd k2 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd k3 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd k4 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(x.rows());
  double t = t0;
  double h = rk_h;
  bool end_flag = false;

  while (!end_flag) {
    if (t + h > t_end) {
      h = t_end - t;
      end_flag = true;
    }
    this->operator()(x, k1, t);
    temp = x + 0.5 * h * k1;
    t += h / 2.0;

    this->operator()(temp, k2, t);
    temp = x + 0.5 * h * k2;

    this->operator()(temp, k3, t);
    temp = x + h * k3;
    t += h / 2.0;

    this->operator()(temp, k4, t);
    x += (h / 6.0) * (k1 + 2.0 * (k2 + k3) + k4);
  }
  // store final state (k1 is put for dummy)
  this->operator()(x, k1, t);
}

void dynamical_system::operator()(const Eigen::VectorXd &x,
                                  Eigen::VectorXd &dxdt, const double t) {
  // store system function f and its derivative df/d*
  sys_func(x, t);
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