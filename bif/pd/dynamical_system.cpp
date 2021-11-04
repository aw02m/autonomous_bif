#include "dynamical_system.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["x0"].size();
  udim = xdim - 1;

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

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["x0"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;
  this->u0 = h(this->x0);

  tauk = Eigen::VectorXd::Zero(period);
  tauk[0] = json["tau"];

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

  dTdx = Eigen::MatrixXd::Zero(xdim, xdim);
  mu = json["sigma"];
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

Eigen::VectorXd dynamical_system::h(const Eigen::VectorXd &x) {
  Eigen::VectorXd ret(udim);

  for (int i = 0; i < xdim; i++) {
    if (i != p_index) {
      if (i < p_index) {
        ret(i) = x(i);
      } else {
        ret(i - 1) = x(i);
      }
    }
  }

  return ret;
}

Eigen::VectorXd dynamical_system::h_inv(const Eigen::VectorXd &u) {
  Eigen::VectorXd ret(xdim);

  for (int i = 0; i < xdim; i++) {
    if (i != p_index) {
      if (i < p_index) {
        ret(i) = u(i);
      } else {
        ret(i) = u(i - 1);
      }
    } else {
      ret(i) = p_place;
    }
  }

  return ret;
}

void dynamical_system::store_states(const Eigen::VectorXd &v) {
  Eigen::VectorXd u = v(Eigen::seqN(0, udim));
  double tau = v(udim);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd O = Eigen::VectorXd::Zero(size_dphidlambda + size_dphidxdx +
                                            size_dphidxdlambda);
  Eigen::VectorXd state(xdim + size_dphidx + size_dphidlambda + size_dphidxdx +
                        size_dphidxdlambda);

  u0 = u;
  tauk(0) = tau;
  xk[0] = h_inv(u);
  p(var_param) = v(xdim);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I, O;
  integrate(0, state, tauk(0));

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

  dTldu = dhdx * dphidx * dh_invdu;
  dTldtau = dhdx * f;
  dTldlambda = dhdx * dphidlambda;
  dqdu = dqdx * dphidx * dh_invdu;
  dqdtau = dqdx * f;
  dqdlambda = dqdx * dphidlambda;
  std::vector<Eigen::MatrixXd> dTdxdu_buf =
      std::vector<Eigen::MatrixXd>(udim, Eigen::MatrixXd::Zero(xdim, xdim));
  for (int i = 0; i < xdim; i++) {
    if (i != p_index) {
      if (i < p_index) {
        dTdxdu_buf[i] = dhdx * dphidxdx[i] * dh_invdu;
      } else {
        dTdxdu_buf[i - 1] = dhdx * dphidxdx[i] * dh_invdu;
      }
    }
  }
  dTdxdu = dTdxdu_buf;
  dTdxdtau = dhdx * dphidxdtau * dh_invdu;
  dTdxdlambda = dhdx * dphidxdlambda * dh_invdu;

  dTdx = dphidx;
  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(dTdx).eigenvalues();
  chara_poly = dTdx - mu * Eigen::MatrixXd::Identity(xdim, xdim);
}

void dynamical_system::store_constant_state() {
  dqdx = Eigen::MatrixXd::Zero(1, xdim);
  dqdx(0, p_index) = 1.0;

  dhdx = Eigen::MatrixXd::Identity(xdim, xdim);
  unsigned int rowToRemove = p_index;
  unsigned int numRows = dhdx.rows() - 1;
  unsigned int numCols = dhdx.cols();
  if (rowToRemove < numRows)
    dhdx.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        dhdx.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
  dhdx.conservativeResize(numRows, numCols);

  dh_invdu = Eigen::MatrixXd::Identity(xdim, xdim);
  unsigned int colToRemove = p_index;
  numRows = dh_invdu.rows();
  numCols = dh_invdu.cols() - 1;
  if (colToRemove < numCols)
    dh_invdu.block(0, colToRemove, numRows, numCols - colToRemove) =
        dh_invdu.block(0, colToRemove + 1, numRows, numCols - colToRemove);
  dh_invdu.conservativeResize(numRows, numCols);
}

Eigen::VectorXd dynamical_system::newton_F() {
  Eigen::VectorXd ret(udim + period + 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);

  ret(Eigen::seqN(0, udim)) = h(xk[period]) - u0;
  ret(udim) = q(xk[period]);
  ret(udim + 1) = chara_poly.determinant();

  return ret;
}

Eigen::MatrixXd dynamical_system::newton_J() {
  Eigen::MatrixXd ret(udim + period + 1, udim + period + 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(udim, udim);
  Eigen::MatrixXd dchidu = Eigen::MatrixXd::Zero(1, udim);
  for (int i = 0; i < udim; i++) {
    dchidu(0, i) = det_derivative(chara_poly, dTdxdu[i]);
  }
  double dchidtau = det_derivative(chara_poly, dTdxdtau);
  double dchidlambda = det_derivative(chara_poly, dTdxdlambda);

  ret(Eigen::seqN(0, udim), Eigen::seqN(0, udim)) = dTldu - I;
  ret(Eigen::seqN(0, udim), udim) = dTldtau;
  ret(Eigen::seqN(0, udim), udim + 1) = dTldlambda;
  ret(udim, Eigen::seqN(0, udim)) = dqdu;
  ret(udim, udim) = dqdtau(0, 0);
  ret(udim, udim + 1) = dqdlambda(0, 0);
  ret(udim + 1, Eigen::seqN(0, udim)) = dchidu;
  ret(udim + 1, udim) = dchidtau;
  ret(udim + 1, udim + 1) = dchidlambda;

  return ret;
}

double dynamical_system::det_derivative(const Eigen::MatrixXd &A,
                                        const Eigen::MatrixXd &dA) {
  Eigen::MatrixXd temp(xdim, xdim);
  double ret = 0;

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
  store_states(v);

  return std::make_tuple(newton_F(), newton_J());
}