#include "newton.hpp"
#include "dynamical_system.hpp"
#include <fstream>

void newton(dynamical_system &ds) {
  Eigen::VectorXd vp(ds.xdim + 1);
  vp << ds.u0, ds.tauk, ds.p(ds.var_param);
  Eigen::VectorXd vn(ds.xdim + 1);
  Eigen::VectorXd F(ds.xdim + 1);
  Eigen::MatrixXd J(ds.xdim + 1, ds.xdim + 1);
  Eigen::VectorXcd eigvals;
  double norm;
  Eigen::IOFormat Comma(8, 0, ", ", "\n", "[", "]");
  std::vector<Eigen::VectorXd> bifset(ds.inc_iter,
                                      Eigen::VectorXd::Zero(ds.p.size()));
  bool exit_flag = false;

  for (int p = 0; p < ds.inc_iter; p++) {
    if (exit_flag) {
      bifset.erase(bifset.begin() + p, bifset.end());
      ds.inc_iter = p - 1;
      break;
    }
    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < ds.max_iter; i++) {
      auto FJ = ds.newton_FJ(vp);
      F = std::get<0>(FJ);
      J = std::get<1>(FJ);
      vn = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(J).solve(-F) + vp;
      norm = (vn - vp).norm();
      // debug(F);
      // debug(J);
      // debug(vp);
      // debug(vn);
      // debug(norm);
      if (norm < ds.eps) {
        auto end = std::chrono::system_clock::now();
        auto dur = end - start;
        auto msec =
            std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        std::cout << "**************************************************"
                  << std::endl;
        std::cout << p << " : converged (iter = " << i + 1 << ", ";
        std::cout << "time = " << msec << "[msec])" << std::endl;
        std::cout << "params : " << ds.p.transpose().format(Comma) << std::endl;
        std::cout << "x0     : "
                  << ds.h_inv(vn(Eigen::seqN(0, ds.udim))).transpose().format(Comma)
                  << std::endl;
        std::cout << "tau    : " << vn(ds.udim) << std::endl;
        std::cout << "(Re(μ), Im(μ)), abs(μ), arg(μ) :" << std::endl;
        for (int k = 0; k < ds.xdim; k++) {
          std::cout << ds.eigvals(k) << ", ";
          std::cout << std::abs(ds.eigvals(k)) << ", ";
          std::cout << std::arg(ds.eigvals(k)) * (180 / EIGEN_PI) << std::endl;
        }
        std::cout << "**************************************************"
                  << std::endl;
        ds.p(ds.var_param) = vn(ds.xdim);
        bifset[p] = ds.p;
        vp = vn;
        break;
      } else if (norm >= ds.explode) {
        std::cerr << "explode (iter = " << i + 1 << ")" << std::endl;
        exit_flag = true;
        break;
      }
      if (i >= ds.max_iter - 1) {
        std::cerr << "iter over (iter = " << i + 1 << ")" << std::endl;
        exit_flag = true;
        break;
      }
      vp = vn;
    }
    ds.p[ds.inc_param] += ds.delta_inc;
  }

  std::ofstream f;
  f.open("out", std::ios::out);
  for (int i = 0; i < ds.inc_iter; i++) {
    f << bifset[i].transpose() << std::endl;
  }
  f.close();
}