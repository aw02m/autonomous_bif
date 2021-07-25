#include "newton.hpp"

void newton(dynamical_system &ds) {
  Eigen::VectorXd vp(ds.xdim);
  vp << ds.u0, ds.tauk;
  Eigen::VectorXd vn(ds.xdim);
  Eigen::VectorXd F(ds.xdim);
  Eigen::MatrixXd J(ds.xdim, ds.xdim);
  Eigen::VectorXcd eigvals(ds.xdim);
  double norm;

  store_constant_state(ds);

  for (int p = 0; p < ds.inc_iter; p++) {
    for (int i = 0; i < ds.max_iter; i++) {
      store_state(vp, ds);

      F = func_newton(ds);
      J = jac_newton(ds);
      vn = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(J).solve(-F) + vp;

      norm = (vn - vp).norm();
      if (norm < ds.eps) {
        std::cout << "******************************************" << std::endl;
        std::cout << p << " : converged (iter = " << i + 1 << ")" << std::endl;
        std::cout << "params : " << ds.params.transpose() << std::endl;
        std::cout << "states : " << vn.transpose() << std::endl;
        eigvals = Eigen::EigenSolver<Eigen::MatrixXd>(ds.dphidx[0]).eigenvalues();
        std::cout << "eigen value of dTldu (Re(μ), Im(μ)) :" << std::endl
                  << eigvals << std::endl;
        std::cout << "abs = " << std::abs(eigvals(0))
                  << ",  arg = " << std::arg(eigvals(0)) << std::endl;
        std::cout << "******************************************" << std::endl;
        vp = vn;
        break;
      } else if (norm >= ds.explode) {
        std::cerr << "explode (iter = " << i + 1 << ")" << std::endl;
        exit(1);
      }

      if (i >= ds.max_iter - 1) {
        std::cerr << "iter over (iter = " << i + 1 << ")" << std::endl;
        exit(1);
      }

      vp = vn;
    }
    ds.params[ds.inc_param] += ds.delta_inc;
  }
}