#include "eigensolver.hpp"
#include "dynamical_system.hpp"

Eigen::VectorXcd eigenvalues(const dynamical_system &ds) {
  Eigen::VectorXcd eigvals(ds.udim);
  // eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(ds.dTldu).eigenvalues();
  eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(ds.dphidx[0]).eigenvalues();

  return eigvals;
}

// take eigenvalue which is nearest to e^j\theta
Eigen::dcomplex bifeigvals(const dynamical_system &ds) {
  unsigned int target_index = 0;
  double delta = std::abs(ds.eigvals(0)) - 1.0;
  double delta_buf = 0;

  for (int i = 1; i < ds.xdim; i++) {
    delta_buf = std::abs(ds.eigvals(i)) - 1.0;
    if (delta_buf < delta) {
      target_index = i;
    }
  }

  return ds.eigvals(target_index);
}