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
  rkf_first_h = json["rkf_first_h"];
  rkf_max_h = json["rkf_max_h"];
  rkf_min_h = json["rkf_min_h"];
  rkf_tol = json["rkf_tol"];
  rkf_false_iter = json["rkf_false_iter"];

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

  integrate(state, 0, tau);

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
  Eigen::VectorXd O = Eigen::VectorXd::Zero(size_dphidlambda);
  Eigen::VectorXd state(xdim + size_dphidx + size_dphidlambda);
  Eigen::VectorXd state_init(xdim + size_dphidx + size_dphidlambda);
  Eigen::VectorXd h = Eigen::VectorXd::Zero(xdim);
  Eigen::MatrixXd Ah;

  xk[0] = v(Eigen::seqN(0, xdim));
  tau = v(xdim);
  p(var_param) = v(xdim + 1);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I, O;
  state_init = state;

  integrate(state, 0, tau);

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
    integrate(state, 0, tau);
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
  integrate(state, 0, tau);
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

void dynamical_system::store_states_fix(const Eigen::VectorXd &v) {
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(xdim, xdim);
  Eigen::VectorXd state(xdim + size_dphidx);

  xk[0] = v(Eigen::seqN(0, xdim));
  tau = v(xdim);

  I.resize(I.cols() * I.rows(), 1);
  state << xk[0], I;

  integrate(state, 0, tau);

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