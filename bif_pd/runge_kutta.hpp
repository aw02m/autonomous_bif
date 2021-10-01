#ifndef RUNGE_KUTTA_HPP_
#define RUNGE_KUTTA_HPP_

#include "sys_common.hpp"

class dynamical_system;

Eigen::VectorXd
integrate(const std::function<Eigen::VectorXd(double, const Eigen::VectorXd &,
                                        const dynamical_system &)>
              func,
          double t0, const Eigen::VectorXd &x, double tau,
          const dynamical_system &ds);

Eigen::VectorXd
integrate_rk45(const std::function<Eigen::VectorXd(double, const Eigen::VectorXd &,
                                        const dynamical_system &)>
              func,
          double t0, const Eigen::VectorXd &x, double tau,
          const dynamical_system &ds);

#endif