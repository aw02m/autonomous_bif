#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

// #include <QString>
#include "qcustomplot.h"
#include <QVector>
#include <eigen3/Eigen/Dense>
#include <vector>

class dynamical_system {
public:
  dynamical_system(const std::string &json_location);

  Eigen::VectorXd func(double t, const Eigen::VectorXd &x);
  double q(const Eigen::VectorXd &x);
  // Eigen::MatrixXd dqdx(const Eigen::VectorXd &x);

  unsigned int xdim;
  // unsigned int p_index;
  // double p_place;
  std::vector<double> q_coef;
  double direction;
  unsigned int period;

  // graphic & rk params
  double tick;
  std::vector<unsigned int> axis;
  std::vector<double> xrange;
  std::vector<double> yrange;
  unsigned int max_plot;
  unsigned int max_poincare_plot;
  bool use_classic_rk;
  double rkf_first_h;
  double rkf_h_max;
  double rkf_h_min;
  double rkf_tol;
  unsigned int rkf_false_iter;
  double eps;
  double explode;
  unsigned int max_iter;
  double delta_inc;

  Eigen::VectorXd x0;
  Eigen::VectorXd p;
  double tau = 0;

  Eigen::VectorXd last_state;
  QVector<QCPCurveData> QCPCsol;
  QVector<QCPGraphData> QCPGpoincare;
  Eigen::MatrixXd dqdx;
  bool hit_section = false;
  bool is_hit_section();

  void integrate(double t0, const Eigen::VectorXd &x, double t_end);
  void integrate_rk45(double t0, const Eigen::VectorXd &x, double t_end);
};

#endif