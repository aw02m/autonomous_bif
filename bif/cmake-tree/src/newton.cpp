#include "newton.hpp"
#include "dynamical_system.hpp"
#include <fstream>
#include <iomanip>

void newton(dynamical_system &ds) {
  unsigned int target_dim;
  switch (ds.mode) {
  case 0:
    target_dim = ds.xdim + 1;
    break;
  case 1:
  case 2:
  case 3:
    target_dim = ds.xdim + 2;
    break;
  case 4:
    target_dim = ds.xdim;
    break;
  case 5:
  case 6:
    target_dim = ds.xdim + 1;
    break;
  }
  Eigen::VectorXd vp(target_dim);
  Eigen::VectorXd vn(target_dim);
  Eigen::VectorXd v_last_succeed(target_dim);
  Eigen::VectorXd F(target_dim);
  Eigen::MatrixXd J(target_dim, target_dim);

  Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]");
  std::cout << std::fixed << std::setprecision(16);

  std::ofstream f;
  Eigen::IOFormat Out(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "",
                      "");
  if (ds.mode != 0 && ds.mode != 4) {
    f.open(ds.out_path, std::ios::out);
    f << std::fixed;
  }

  bool exit_flag = false;

  double norm = 1;
  switch (ds.mode) {
  case 0:
    vp(Eigen::seqN(0, ds.xdim)) = ds.x0;
    vp(ds.xdim) = ds.tau;
    break;
  case 1:
  case 2:
  case 3:
    vp(Eigen::seqN(0, ds.xdim)) = ds.x0;
    vp(ds.xdim) = ds.tau;
    vp(ds.xdim + 1) = ds.p(ds.var_param);
    break;
  case 4:
    vp(Eigen::seqN(0, ds.xdim)) = ds.x0;
    break;
  case 5:
  case 6:
    vp(Eigen::seqN(0, ds.xdim)) = ds.x0;
    vp(ds.xdim) = ds.p(ds.var_param);
    break;
  }

  for (int p = 0; p < ds.inc_iter; p++) {
    if (exit_flag) {
      ds.inc_iter = p - 1;
      break;
    }
    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < ds.max_iter; i++) {
      auto FJ = ds.newton_FJ(vp);
      F = std::get<0>(FJ);
      J = std::get<1>(FJ);
      vn = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(J).solve(-F) + vp;
      // std::cout << "error  : " << F.transpose().format(Comma) << std::endl;
      norm = F.norm();
      if (check_norm(F, ds)) {
        // debug_exit(F);
        auto end = std::chrono::system_clock::now();
        auto dur = end - start;
        auto msec =
            std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

        std::cout << std::setprecision(16);
        std::cout << p << " : converged (iter = " << i + 1 << ", ";
        std::cout << "error-norm = " << norm << ", ";
        std::cout << "time = " << msec << "[msec])" << std::endl;
        // std::cout << "error  : " << F.transpose().format(Comma) << std::endl;
        std::cout << "x0     : "
                  << vn(Eigen::seqN(0, ds.xdim)).transpose().format(Comma)
                  << std::endl;
        std::cout << "tau    : " << vn(ds.xdim) << std::endl;
        std::cout << "params : " << ds.p.transpose().format(Comma) << std::endl;
        // std::cout << "theta  : " << ds.theta << std::endl;
        std::cout << "(Re(μ), Im(μ)), abs(μ), arg(μ) :" << std::endl;
        for (int k = 0; k < ds.xdim; k++) {
          std::cout << std::setprecision(8);
          std::cout << ds.eigvals(k) << ", ";
          std::cout << std::abs(ds.eigvals(k)) << ", ";
          std::cout << std::arg(ds.eigvals(k)) << std::endl;
        }
        std::cout << "**************************************************"
                  << std::endl;
        f << ds.p.transpose().format(Out) << std::endl;
        v_last_succeed = vn;
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

  // set last state
  switch (ds.mode) {
  case 0:
    ds.x0 = v_last_succeed(Eigen::seqN(0, ds.xdim));
    ds.tau = v_last_succeed(ds.xdim);
    break;
  case 1:
  case 2:
  case 3:
    ds.x0 = v_last_succeed(Eigen::seqN(0, ds.xdim));
    ds.tau = v_last_succeed(ds.xdim);
    ds.p(ds.var_param) = v_last_succeed(ds.xdim + 1);
    break;
  case 4:
    ds.x0 = v_last_succeed(Eigen::seqN(0, ds.xdim));
    break;
  case 5:
  case 6:
    ds.x0 = v_last_succeed(Eigen::seqN(0, ds.xdim));
    ds.p(ds.var_param) = v_last_succeed(ds.xdim);
    break;
  }

  if (ds.mode != 0 && ds.mode != 4)
    f.close();
}

bool check_norm(Eigen::VectorXd &F, dynamical_system &ds) {
  bool ret = true;
  unsigned int dim = F.rows();
  for (int i = 0; i < dim; i++) {
    if (std::abs(F(i)) > ds.eps)
      ret = false;
  }
  return ret;
}