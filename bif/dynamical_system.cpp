#include "dynamical_system.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["x0"].size();

  p_index = json["p_index"];
  p_place = json["p_place"];

  use_classic_rk = json["use_classic_rk"];
  rk_div = json["rk_div"];
  rkf_first_h = json["rkf_first_h"];
  rkf_h_max = json["rkf_h_max"];
  rkf_h_min = json["rkf_h_min"];
  rkf_tol = json["rkf_tol"];
  rkf_false_iter = json["rkf_false_iter"];

  period = json["period"];
  inc_param = json["inc_param"];
  var_param = json["var_param"];
  delta_inc = json["delta_inc"];
  inc_iter = json["inc_iter"];
  max_iter = json["max_iter"];
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
  theta = json["theta"];

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

  switch (mode) {
  case 1:
    mu = 1;
    break;
  case 2:
    mu = -1;
    break;
  case 3:
    mu = Eigen::dcomplex(std::cos(theta), std::sin(theta));
    break;
  }
}

void dynamical_system::integrate(double t_0, Eigen::VectorXd &x, double t_end) {
  Eigen::VectorXd k1 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd k2 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd k3 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd k4 = Eigen::VectorXd::Zero(x.rows());
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(x.rows());
  double t = t_0;
  double h = t_end / rk_div;

  for (int i = 0; i < rk_div; i++) {
    function(t, x, k1);
    temp = x + 0.5 * h * k1;
    t += h / 2.0;

    function(t, temp, k2);
    temp = x + 0.5 * h * k2;

    function(t, temp, k3);
    temp = x + h * k3;
    t += h / 2.0;

    function(t, temp, k4);
    x += (h / 6.0) * (k1 + 2.0 * (k2 + k3) + k4);
  }
  // store final state (k1 is put for dummy)
  function(t, x, k1);
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

  switch (mode) {
  case 3:
    theta = v(xdim + 2);
  }

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I, O;
  integrate(0, state, tau);

  // store derivatives of final state
  Eigen::VectorXd dummy(xdim + size_dphidx + size_dphidlambda + size_dphidxdx +
                        size_dphidxdlambda);
  function(0, state, dummy);

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
  dqdx = dqdx * dphidx;
  dqdtau = dqdx * f;
  dqdlambda = dqdx * dphidlambda;
  dTdxdx = dphidxdx;
  dTdxdtau = dphidxdtau;
  dTdxdlambda = dphidxdlambda;

  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTdx).eigenvalues();

  if (mode == 3)
    mu = Eigen::dcomplex(std::cos(theta), std::sin(theta));

  I = Eigen::MatrixXd::Identity(xdim, xdim);
  chara_poly = dTdx - mu * Eigen::MatrixXd::Identity(xdim, xdim);
}

void dynamical_system::store_constant_state() {
  dqdx = Eigen::MatrixXd::Zero(1, xdim);
  dqdx(0, p_index) = 1.0;
}

Eigen::VectorXd dynamical_system::newton_F_PD() {
  Eigen::VectorXd ret(xdim + 2);
  double chi = chara_poly.determinant().real();

  ret(Eigen::seqN(0, xdim)) = xk[1] - xk[0];
  ret(xdim) = q(xk[1]);
  ret(xdim + 1) = chi;

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_PD() {
  Eigen::MatrixXd ret(xdim + 2, xdim + 2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::MatrixXd dchidx = Eigen::MatrixXd::Zero(1, xdim);
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dTdxdx[i]).real();
  }
  double dchidtau = det_derivative(chara_poly, dTdxdtau).real();
  double dchidlambda = det_derivative(chara_poly, dTdxdlambda).real();

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTdx - I;
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(Eigen::seqN(0, xdim), xdim + 1) = dTdlambda;
  ret(xdim, Eigen::seqN(0, xdim)) = dqdx;
  ret(xdim, xdim) = dqdtau(0, 0);
  ret(xdim, xdim + 1) = dqdlambda(0, 0);
  ret(xdim + 1, Eigen::seqN(0, xdim)) = dchidx;
  ret(xdim + 1, xdim) = dchidtau;
  ret(xdim + 1, xdim + 1) = dchidlambda;

  return ret;
}

Eigen::VectorXd dynamical_system::newton_F_NS() {
  Eigen::VectorXd ret(xdim + 3);
  Eigen::dcomplex chi = chara_poly.determinant().real();

  ret(Eigen::seqN(0, xdim)) = xk[1] - xk[0];
  ret(xdim) = q(xk[1]);
  ret(xdim + 1) = chi.real();
  ret(xdim + 1) = chi.imag();

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_NS() {
  Eigen::MatrixXd ret(xdim + 3, xdim + 3);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::MatrixXcd dchidx = Eigen::MatrixXcd::Zero(1, xdim);
  for (int i = 0; i < xdim; i++) {
    dchidx(0, i) = det_derivative(chara_poly, dTdxdx[i]);
  }
  Eigen::dcomplex dchidtau = det_derivative(chara_poly, dTdxdtau);
  Eigen::dcomplex dchidlambda = det_derivative(chara_poly, dTdxdlambda);
  Eigen::MatrixXcd dpolydtheta =
      Eigen::dcomplex(std::sin(theta), -std::cos(theta)) * I;
  Eigen::dcomplex dchidtheta = det_derivative(chara_poly, dpolydtheta);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTdx - I;
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(Eigen::seqN(0, xdim), xdim + 1) = dTdlambda;
  ret(Eigen::seqN(0, xdim), xdim + 2) = Eigen::VectorXd::Zero(xdim);
  ret(xdim, Eigen::seqN(0, xdim)) = dqdx;
  ret(xdim, xdim) = dqdtau(0, 0);
  ret(xdim, xdim + 1) = dqdlambda(0, 0);
  ret(xdim, xdim + 2) = 0;
  ret(xdim + 1, Eigen::seqN(0, xdim)) = dchidx.real();
  ret(xdim + 1, xdim) = dchidtau.real();
  ret(xdim + 1, xdim + 1) = dchidlambda.real();
  ret(xdim + 1, xdim + 2) = dchidtheta.real();
  ret(xdim + 2, Eigen::seqN(0, xdim)) = dchidx.imag();
  ret(xdim + 2, xdim) = dchidtau.imag();
  ret(xdim + 2, xdim + 1) = dchidlambda.imag();
  ret(xdim + 2, xdim + 2) = dchidtheta.imag();

  return ret;
}

void dynamical_system::store_states_fix(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd state(xdim + size_dphidx);

  xk[0] = v(Eigen::seqN(0, xdim));
  tau = v(xdim);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I;
  integrate(0, state, tau);

  // store derivatives of final state
  Eigen::VectorXd dummy(xdim + size_dphidx + size_dphidlambda + size_dphidxdx +
                        size_dphidxdlambda);
  function(0, state, dummy);

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
  dqdx = dqdx * dphidx;
  dqdtau = dqdx * f;

  // Find the argument <theta> of the characteristic constant
  // whose absolute value is closest to 1.
  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTdx).eigenvalues();
  unsigned int target_index = 0;
  double delta = std::abs(eigvals(0)) - 1.0;
  double delta_buf = 0;
  for (int i = 1; i < xdim; i++) {
    delta_buf = std::abs(eigvals(i)) - 1.0;
    if (delta_buf < delta) {
      target_index = i;
    }
  }
  theta = std::arg(eigvals(target_index));
}

Eigen::VectorXd dynamical_system::newton_F_fix() {
  Eigen::VectorXd ret(xdim + 1);

  ret(Eigen::seqN(0, xdim)) = xk[1] - xk[0];
  ret(xdim) = q(xk[1]);

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J_fix() {
  Eigen::MatrixXd ret(xdim + 1, xdim + 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);

  ret(Eigen::seqN(0, xdim), Eigen::seqN(0, xdim)) = dTdx - I;
  ret(Eigen::seqN(0, xdim), xdim) = dTdtau;
  ret(xdim, Eigen::seqN(0, xdim)) = dqdx;
  ret(xdim, xdim) = dqdtau(0, 0);

  return ret;
}

Eigen::dcomplex dynamical_system::det_derivative(const Eigen::MatrixXcd &A,
                                                 const Eigen::MatrixXcd &dA) {
  Eigen::MatrixXcd temp(xdim, xdim);
  Eigen::dcomplex ret = 0;

  for (int i = 0; i < xdim; i++) {
    temp = A;
    temp.col(i) = dA.col(i);
    ret += temp.determinant();
  }

  return ret;
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
dynamical_system::newton_FJ(const Eigen::VectorXd &v) {
  store_constant_state();
  switch (mode) {
  case 0:
    store_states_fix(v);
    return std::make_tuple(newton_F_fix(), newton_J_fix());
  case 2:
    store_states(v);
    return std::make_tuple(newton_F_PD(), newton_J_PD());
  case 3:
    store_states(v);
    return std::make_tuple(newton_F_NS(), newton_J_NS());
  }
  // exception
  return std::make_tuple(Eigen::VectorXd::Zero(xdim),
                         Eigen::MatrixXd::Zero(xdim, xdim));
}