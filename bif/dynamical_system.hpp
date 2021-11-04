#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include "sys_common.hpp"
#include <nlohmann/json.hpp>

class dynamical_system {
public:
  dynamical_system(nlohmann::json json);
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  newton_FJ(const Eigen::VectorXd &v);

  unsigned int mode;
  std::string out_path;
  std::string json_out_path;

  unsigned int xdim;
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
  double tau;
  Eigen::VectorXd p;
  double theta;

  Eigen::MatrixXcd eigvals;

private:
  std::vector<Eigen::VectorXd> xk;

  Eigen::VectorXd f;
  Eigen::MatrixXd dfdx;
  Eigen::VectorXd dfdlambda;
  std::vector<Eigen::MatrixXd> dfdxdx;
  Eigen::MatrixXd dfdxdlambda;


  Eigen::MatrixXd dphidx;
  Eigen::VectorXd dphidlambda;
  std::vector<Eigen::MatrixXd> dphidxdx;
  Eigen::MatrixXd dphidxdtau;
  Eigen::MatrixXd dphidxdlambda;

  unsigned int size_dphidx;
  unsigned int size_dphidlambda;
  unsigned int size_dphidxdx;
  unsigned int size_dphidxdlambda;

  Eigen::MatrixXd dTdx;
  Eigen::VectorXd dTdtau;
  Eigen::VectorXd dTdlambda;
  Eigen::MatrixXd dqdx;
  Eigen::MatrixXd dqdtau;
  Eigen::MatrixXd dqdlambda;
  std::vector<Eigen::MatrixXd> dTdxdx;
  Eigen::MatrixXd dTdxdtau;
  Eigen::MatrixXd dTdxdlambda;

  Eigen::dcomplex mu;
  Eigen::MatrixXcd chara_poly;

  void function(double t, const Eigen::VectorXd &x, Eigen::VectorXd &dxdt);
  double q(const Eigen::VectorXd &x);
  void integrate(double t_0, Eigen::VectorXd &x, double t_end);
  Eigen::dcomplex det_derivative(const Eigen::MatrixXcd &A,
                                 const Eigen::MatrixXcd &dA);

  void store_constant_state();
  void store_states(const Eigen::VectorXd &v);
  void store_states_fix(const Eigen::VectorXd &v);
  Eigen::VectorXd newton_F_fix();
  Eigen::MatrixXd newton_J_fix();
  Eigen::VectorXd newton_F_G();
  Eigen::MatrixXd newton_J_G();
  Eigen::VectorXd newton_F_PD();
  Eigen::MatrixXd newton_J_PD();
  Eigen::VectorXd newton_F_NS();
  Eigen::MatrixXd newton_J_NS();
};

#endif