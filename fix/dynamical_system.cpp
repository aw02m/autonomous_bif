#include "dynamical_system.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["fixed"].size();
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
  dif_strip = json["dif_strip"];
  eps = json["eps"];
  explode = json["explode"];

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["fixed"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;
  this->u0 = h(this->x0, *this);

  tauk = Eigen::VectorXd::Zero(period);
  tauk[0] = json["tau"];

  std::vector<double> params_arr = json["params"];
  Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
  this->params = params;

  xk = std::vector<Eigen::VectorXd>(period + 1, Eigen::VectorXd::Zero(xdim));
  fk = std::vector<Eigen::VectorXd>(period, Eigen::VectorXd::Zero(xdim));
  dfdx =
      std::vector<Eigen::MatrixXd>(period, Eigen::MatrixXd::Zero(xdim, xdim));
  dphidx =
      std::vector<Eigen::MatrixXd>(period, Eigen::MatrixXd::Zero(xdim, xdim));
}