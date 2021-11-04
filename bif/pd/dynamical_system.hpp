#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include "sys_common.hpp"
#include <nlohmann/json.hpp>

class dynamical_system {
public:
  dynamical_system(nlohmann::json json);
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  newton_FJ(const Eigen::VectorXd &v);

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
  double eps;
  double explode;

  Eigen::VectorXd x0;
  Eigen::VectorXd u0;
  Eigen::VectorXd tauk;
  Eigen::VectorXd p;

  Eigen::MatrixXcd eigvals;
  
  Eigen::VectorXd h_inv(const Eigen::VectorXd &u);

private:
  std::vector<Eigen::VectorXd> xk;

  Eigen::VectorXd f;
  Eigen::MatrixXd dfdx;
  Eigen::VectorXd dfdlambda;
  std::vector<Eigen::MatrixXd> dfdxdx;
  Eigen::MatrixXd dfdxdlambda;  

  Eigen::MatrixXd dqdx;
  Eigen::MatrixXd dhdx;
  Eigen::MatrixXd dh_invdu;

  Eigen::MatrixXd dphidx;
  Eigen::VectorXd dphidlambda;
  std::vector<Eigen::MatrixXd> dphidxdx;
  Eigen::MatrixXd dphidxdtau;
  Eigen::MatrixXd dphidxdlambda;

  unsigned int size_dphidx;
  unsigned int size_dphidlambda;
  unsigned int size_dphidxdx;
  unsigned int size_dphidxdlambda;

  Eigen::MatrixXd chara_poly;

  Eigen::MatrixXd dTldu;
  Eigen::VectorXd dTldtau;
  Eigen::VectorXd dTldlambda;
  Eigen::MatrixXd dqdu;
  Eigen::MatrixXd dqdtau;
  Eigen::MatrixXd dqdlambda;
  std::vector<Eigen::MatrixXd> dTdxdu;
  Eigen::MatrixXd dTdxdtau;
  Eigen::MatrixXd dTdxdlambda;

  Eigen::MatrixXd dTdx;
  double mu;

  void function(double t, const Eigen::VectorXd &x, Eigen::VectorXd &dxdt);
  double q(const Eigen::VectorXd &x);
  Eigen::VectorXd h(const Eigen::VectorXd &x);
  void integrate(double t_0, Eigen::VectorXd &x, double t_end);
  double det_derivative(const Eigen::MatrixXd &A, const Eigen::MatrixXd &dA);

  void store_constant_state();
  void store_states(const Eigen::VectorXd &v);
  Eigen::VectorXd newton_F();
  Eigen::MatrixXd newton_J();
};

#endif