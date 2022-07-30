#include "dynamical_system.hpp"

void dynamical_system::integrate(Eigen::VectorXd &x, double t0, double t_end) {
  if (use_classic_rk) {
    rk45(x, t0, t_end);
  } else {
    rkf45(x, t0, t_end);
  }
}

void dynamical_system::rk45(Eigen::VectorXd &x, double t0, double t_end) {
  unsigned int dim = x.rows();
  Eigen::VectorXd k1 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k2 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k3 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k4 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(dim);
  double t = t0;
  double h = rkf_first_h;
  bool end_flag = false;

  while (!end_flag) {
    if (t + h > t_end) {
      h = t_end - t;
      end_flag = true;
    }
    this->operator()(x, k1, t);
    temp = x + 0.5 * h * k1;
    t += h / 2.0;

    this->operator()(temp, k2, t);
    temp = x + 0.5 * h * k2;

    this->operator()(temp, k3, t);
    temp = x + h * k3;
    t += h / 2.0;

    this->operator()(temp, k4, t);
    x += (h / 6.0) * (k1 + 2.0 * (k2 + k3) + k4);
  }
  // store final state (k1 is put for dummy)
  this->operator()(x, k1, t);
}

void dynamical_system::rkf45(Eigen::VectorXd &x, double t0, double t_end) {
  unsigned int dim = x.rows();
  Eigen::VectorXd k1 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k2 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k3 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k4 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k5 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd k6 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(dim);
  double t = t0;
  double h = rkf_first_h;
  double h_max = rkf_max_h;
  double h_min = rkf_min_h;
  double tol = rkf_tol;
  double R = 0.0;
  double delta = 0.0;
  unsigned int loop_counter = 0;
  bool flag = true;

  while (flag) {
    this->operator()(x, k1, t);
    k1 *= h;
    temp = x + k1 * 0.25;
    this->operator()(temp, k2, t + h * 0.25);
    k2 *= h;
    // k2 = h * this->operator()(t0 + h * 0.25, temp, t);
    temp = x + 0.09375 * k1 + 0.28125 * k2;
    this->operator()(temp, k3, t + 0.375 * h);
    k3 *= h;
    // k3 = h * this->operator()(t0 + 0.375 * h, temp, t);
    temp = x + 0.8793809740555303 * k1 - 3.277196176604461 * k2 +
           3.3208921256258535 * k3;
    this->operator()(temp, k4, t + 0.9230769230769231 * h);
    k4 *= h;
    // k4 = h * this->operator()(t0 + 0.9230769230769231 * h, temp, t);
    temp = x + 2.0324074074074074 * k1 - 8.0 * k2 + 7.173489278752436 * k3 -
           0.20589668615984405 * k4;
    this->operator()(temp, k5, t + h);
    k5 *= h;
    // k5 = h * this->operator()(t0 + h, temp, t);
    temp = x - 0.2962962962962963 * k1 + 2.0 * k2 - 1.3816764132553607 * k3 +
           0.4529727095516569 * k4 - 0.275 * k4;
    this->operator()(temp, k6, t + 0.5 * h);
    // k6 = h * this->operator()(t0 + 0.5 * h, temp, t);

    R = (0.002777777777777778 * k1 - 0.02994152046783626 * k3 -
         0.029199893673577886 * k4 + 0.02 * k5 + 0.03636363636363636 * k6)
            .norm() /
        h;

    if (R <= tol) {
      loop_counter = 0;
      t += h;
      x += 0.11574074074074074 * k1 + 0.5489278752436647 * k3 +
           0.5353313840155945 * k4 - 0.2 * k5;
    }

    if (loop_counter++ > rkf_false_iter) {
      std::cerr << "RKF5(6) got stucked. Exitting." << std::endl;
      std::cerr << "Check your system function and derivative." << std::endl;
      exit(1);
    }

    delta = std::pow((tol / (2.0 * R)), 0.25);
    if (delta <= 0.1) {
      h *= 0.1;
    } else if (delta >= 4) {
      h *= 4.0;
    } else {
      h *= delta;
    }

    if (h >= h_max) {
      h = h_max;
    }
    if (h <= h_min) {
      h = h_min;
    }

    if (t >= t_end) {
      flag = false;
    } else if (t + h > t_end) {
      h = t_end - t;
      if (h < h_min) {
        flag = false;
      }
    }
  }
}