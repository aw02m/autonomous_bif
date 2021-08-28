#ifndef DS_DERIVATIVES_HPP_
#define DS_DERIVATIVES_HPP_

#include "sys_common.hpp"

class dynamical_system;

Eigen::VectorXd f(double t, const Eigen::VectorXd &x,
                  const dynamical_system &ds);
Eigen::MatrixXd dfdx(const Eigen::VectorXd &x, const dynamical_system &ds);
Eigen::VectorXd dfdlambda(const Eigen::VectorXd &x, const dynamical_system &ds);
std::vector<Eigen::MatrixXd> dfdxdx(const Eigen::VectorXd &x, const dynamical_system &ds);
Eigen::MatrixXd dfdxdlambda(const Eigen::VectorXd &x,
                            const dynamical_system &ds);

#endif