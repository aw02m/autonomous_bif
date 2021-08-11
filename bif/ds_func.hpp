#ifndef DS_FUNC_HPP_
#define DS_FUNC_HPP_

#include "sys_common.hpp"

class dynamical_system;

void store_constant_state(dynamical_system &ds);
void store_state(const Eigen::VectorXd &vp, dynamical_system &ds);

Eigen::VectorXd func_newton(const dynamical_system &ds);
Eigen::MatrixXd jac_newton(const dynamical_system &ds);

Eigen::VectorXd T(const Eigen::VectorXd &x, double tau,
                  const dynamical_system &ds);

double q(const Eigen::VectorXd &x, const dynamical_system &ds);
Eigen::VectorXd h(const Eigen::VectorXd &x, const dynamical_system &ds);
Eigen::VectorXd h_inv(const Eigen::VectorXd &u, const dynamical_system &ds);
Eigen::MatrixXd dqdx(const dynamical_system &ds);
Eigen::MatrixXd dhdx(const dynamical_system &ds);
Eigen::MatrixXd dh_invdu(const dynamical_system &ds);

Eigen::VectorXd variational_eq(double t, const Eigen::VectorXd &x,
                               const dynamical_system &ds);
Eigen::VectorXd variational_init1(double t, const Eigen::MatrixXd &dphidx,
                                  const Eigen::MatrixXd &dfdx,
                                  const dynamical_system &ds);
Eigen::VectorXd variational_para1(double t, const Eigen::VectorXd &dphidlambda,
                                  const Eigen::MatrixXd &dfdx,
                                  const dynamical_system &ds);

Eigen::MatrixXd dTldu(const dynamical_system &ds);
Eigen::VectorXd dTldtau(const dynamical_system &ds);
Eigen::VectorXd dTldlambda(const dynamical_system &ds);
Eigen::MatrixXd dqdu(const dynamical_system &ds);
Eigen::MatrixXd dqdtau(const dynamical_system &ds);
Eigen::VectorXd dqdlambda(const dynamical_system &ds);
Eigen::MatrixXd dchidu(const dynamical_system &ds);
double dchidtau(const dynamical_system &ds);
double dchidlambda(const dynamical_system &ds);

double det_derivative(const Eigen::MatrixXd &A, const Eigen::MatrixXd &dA,
                      const dynamical_system &ds);

void removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove);
void removeCol(Eigen::MatrixXd &matrix, unsigned int colToRemove);

#endif