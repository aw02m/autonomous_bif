#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include "sys_common.hpp"

class dynamical_system {
public:
  dynamical_system(nlohmann::json json);

  unsigned int xdim;
  unsigned int udim;
  unsigned int period;
  unsigned int p_index;
  double p_place;

  // RK param
  bool use_classic_rk;
  unsigned int rk_div;
  double rkf_first_h;
  double rkf_h_max;
  double rkf_h_min;
  double rkf_tol;
  unsigned int rkf_false_iter;

  // Newton param
  unsigned int inc_param;
  unsigned int var_param;
  double delta_inc;
  unsigned int inc_iter;
  unsigned int max_iter;
  double dif_strip;
  double eps;
  double explode;

  Eigen::VectorXd x0;
  Eigen::VectorXd u0;
  Eigen::VectorXd tauk;
  Eigen::VectorXd params;

  std::vector<Eigen::VectorXd> xk;
  std::vector<Eigen::VectorXd> fk;

  Eigen::MatrixXd dhdx;
  Eigen::MatrixXd dh_invdu;
  Eigen::MatrixXd dqdx;

  std::vector<Eigen::MatrixXd> dphidx;

  Eigen::MatrixXd dTldu;
};

#endif