#include "newton.hpp"
#include "ds_func.hpp"
#include "dynamical_system.hpp"

void newton(dynamical_system &ds) {
  Eigen::VectorXd vp(ds.udim + ds.period + 1 + 1);
  vp << ds.u0, ds.tauk, ds.params(ds.var_param), ds.theta;
  Eigen::VectorXd vn(ds.udim + ds.period + 1 + 1);
  Eigen::VectorXd F(ds.udim + ds.period + 1 + 1);
  Eigen::MatrixXd J(ds.udim + ds.period + 1 + 1, ds.udim + ds.period + 1 + 1);
  Eigen::VectorXcd eigvals;
  double norm;

  store_constant_state(ds);

  for (int p = 0; p < ds.inc_iter; p++) {
    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < ds.max_iter; i++) {
      store_state(vp, ds);

      F = func_newton(ds);
      J = jac_newton(ds);
      vn = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(J).solve(-F) + vp;

      norm = (vn - vp).norm();
      if (norm < ds.eps) {
        auto end = std::chrono::system_clock::now();
        auto dur = end - start;
        auto msec =
            std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        std::cout << "**************************************************"
                  << std::endl;
        std::cout << p << " : converged (iter = " << i + 1 << ", ";
        std::cout << "time = " << msec << "[msec])" << std::endl;
        std::cout << "params : " << ds.params.transpose() << std::endl;
        std::cout << "u0     : " << vn(Eigen::seqN(0, ds.udim)).transpose()
                  << std::endl;
        std::cout << "tau    : " << vn(ds.udim) << std::endl;
        std::cout << "(Re(μ), Im(μ)), abs(μ), arg(μ) :" << std::endl;
        for (int k = 0; k < ds.udim; k++) {
          std::cout << ds.eigvals(k) << ", ";
          std::cout << std::abs(ds.eigvals(k)) << ", ";
          std::cout << std::arg(ds.eigvals(k)) * (180 / EIGEN_PI) << std::endl;
        }
        std::cout << "**************************************************"
                  << std::endl;
        vp = vn;
        ds.params(ds.var_param) = vn(ds.udim + ds.period);
        ds.theta = vn(ds.udim + ds.period + 1);
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
      ds.params(ds.var_param) = vn(ds.udim + ds.period);
      ds.theta = vn(ds.udim + ds.period + 1);
    }
    ds.params[ds.inc_param] += ds.delta_inc;
  }
}