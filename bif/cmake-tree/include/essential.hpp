#pragma once

#include <eigen3/Eigen/Dense>

double det_derivative(const Eigen::MatrixXd &A, const Eigen::MatrixXd &dA,
                      unsigned int dim);
Eigen::MatrixXd bialt_prod(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                           unsigned int xdim, unsigned int bialt_dim);
Eigen::MatrixXd biproduct(const Eigen::MatrixXd &A, unsigned int xdim,
                          unsigned int bialt_dim);
Eigen::MatrixXd biproduct_derivative(const Eigen::MatrixXd &dA,
                                     unsigned int xdim, unsigned int bialt_dim);
Eigen::MatrixXd bialt_prod_square(const Eigen::MatrixXd &A, unsigned int xdim,
                                  unsigned int bialt_dim);
Eigen::MatrixXd bialt_prod_square_derivative(const Eigen::MatrixXd &A,
                                             const Eigen::MatrixXd &dA,
                                             unsigned int xdim,
                                             unsigned int bialt_dim);
unsigned int nearest_idx(const Eigen::VectorXcd &array, double target_val);