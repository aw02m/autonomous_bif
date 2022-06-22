#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include <QVector>

class dynamical_system;

class odesolver : public dynamical_system {
public:
  odesolver(const std::string &json_location);
  void integrate(double t0, const Eigen::VectorXd &x, double t_end);

private:
  unsigned int max_plot;
  unsigned int max_poincare_plot;
  bool use_classic_rk;
  unsigned int rk_div;
  double rkf_first_h;
  double rkf_h_max;
  double rkf_h_min;
  double rkf_tol;
  unsigned int rkf_false_iter;
}

#endif