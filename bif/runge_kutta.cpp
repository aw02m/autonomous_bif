#include "runge_kutta.hpp"

// Runge-Kutta-Fehlberg Method (with automatic step "h" control)
Eigen::VectorXd
integrate(const std::function<Eigen::VectorXd(double, const Eigen::VectorXd &,
                                              const dynamical_system &)>
              func,
          double t0, const Eigen::VectorXd &x, double tau,
          const dynamical_system &ds) {
  Eigen::VectorXd ret = x;
  Eigen::MatrixXd k = Eigen::MatrixXd::Zero(ret.rows(), 6);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(ret.rows());
  double t = t0;
  double h = ds.rkf_first_h;
  double h_max = ds.rkf_h_max;
  double h_min = ds.rkf_h_min;
  double tol = ds.rkf_tol;
  double R = 0.0;
  double delta = 0.0;
  unsigned int counter = 0;
  bool flag = true;

  while (flag) {
    k.col(0) = h * func(t0, ret, ds);
    temp = ret + k.col(0) * 0.25;
    k.col(1) = h * func(t0 + h * 0.25, temp, ds);
    temp = ret + 0.09375 * k.col(0) + 0.28125 * k.col(1);
    k.col(2) = h * func(t0 + 0.375 * h, temp, ds);
    temp = ret + 0.8793809740555303 * k.col(0) - 3.277196176604461 * k.col(1) +
           3.3208921256258535 * k.col(2);
    k.col(3) = h * func(t0 + 0.9230769230769231 * h, temp, ds);
    temp = ret + 2.0324074074074074 * k.col(0) - 8.0 * k.col(1) +
           7.173489278752436 * k.col(2) - 0.20589668615984405 * k.col(3);
    k.col(4) = h * func(t0 + h, temp, ds);
    temp = ret - 0.2962962962962963 * k.col(0) + 2.0 * k.col(1) -
           1.3816764132553607 * k.col(2) + 0.4529727095516569 * k.col(3) -
           0.275 * k.col(3);
    k.col(5) = h * func(t0 + 0.5 * h, temp, ds);

    R = (0.002777777777777778 * k.col(0) - 0.02994152046783626 * k.col(2) -
         0.029199893673577886 * k.col(3) + 0.02 * k.col(4) +
         0.03636363636363636 * k.col(5))
            .norm() /
        h;

    if (R <= tol) {
      counter = 0;
      t += h;
      ret += 0.11574074074074074 * k.col(0) + 0.5489278752436647 * k.col(2) +
             0.5353313840155945 * k.col(3) - 0.2 * k.col(4);
    }
    if (counter++ > ds.rkf_false_iter) {
      std::cerr << "RKF5(6) got stucked. Exitting." << std::endl;
      std::cerr << "Check your system function and derivative." << std::endl;
      exit(1);
    }

    delta = std::pow((tol / (2.0 * R)), 0.25);
    if (delta <= 0.1) {
      h *= 0.1;
    } else if (delta >= 1.5) {
      h *= 2.0;
    } else {
      h *= delta;
    }

    if (h >= h_max) {
      h = h_max;
    }
    if (h <= h_min) {
      h = h_min;
    }

    if (t >= tau) {
      flag = false;
    } else if (t + h > tau) {
      h = tau - t;
      if (h < h_min) {
        flag = false;
      }
    }
  }

  return ret;
}

// Classic Runge-Kutta method
Eigen::VectorXd integrate_rk45(
    const std::function<Eigen::VectorXd(double, const Eigen::VectorXd &,
                                        const dynamical_system &)>
        func,
    double t0, const Eigen::VectorXd &x, double tau,
    const dynamical_system &ds) {
  Eigen::VectorXd ret = x;
  Eigen::MatrixXd k(x.rows(), 4);
  Eigen::VectorXd temp(x.rows());
  double h = tau / ds.rk_div;

  for (int i = 0; i < ds.rk_div; i++) {
    k.col(0) = func(t0, ret, ds);
    temp = ret + h * 0.5 * k.col(0);

    k.col(1) = func(t0 + h * 0.5, temp, ds);
    temp = ret + h * 0.5 * k.col(1);

    k.col(2) = func(t0 + h * 0.5, temp, ds);
    temp = ret + h * k.col(2);

    k.col(3) = func(t0 + h, temp, ds);
    ret += (h / 6.0) * (k.col(0) + 2.0 * k.col(1) + 2.0 * k.col(2) + k.col(3));
  }

  return ret;
}

// // wrapper
// Eigen::VectorXd
// integrate(std::function<Eigen::VectorXd(double, Eigen::VectorXd &,
//                                         dynamical_system &)>
//               func,
//           Eigen::VectorXd &x, double tau, dynamical_system &ds) {
//   Eigen::VectorXd ret = x;
//   double h = tau / ds.rk_div;

//   for (int i = 0; i < ds.rk_div; i++) {
//     ret = runge_kutta(func, 0, ret, h, ds);
//   }

//   return ret;
// }

// // single step of rk45
// Eigen::VectorXd
// runge_kutta(std::function<Eigen::VectorXd(double, Eigen::VectorXd &,
//                                           dynamical_system &)>
//                 func,
//             double t, Eigen::VectorXd &x, double h, dynamical_system &ds) {
//   unsigned int xdim = ds.xdim;
//   Eigen::MatrixXd k(xdim, 4);
//   Eigen::VectorXd temp(xdim);
//   Eigen::VectorXd ret = x;
//   int i;

//   k.col(0) = func(t, x, ds);
//   temp = x + h / 2 * k.col(0);

//   k.col(1) = func(t + h / 2, temp, ds);
//   temp = x + h / 2 * k.col(1);

//   k.col(2) = func(t + h / 2, temp, ds);
//   temp = x + h * k.col(2);

//   k.col(3) = func(t + h, temp, ds);
//   ret += (k.col(0) + 2 * k.col(1) + 2 * k.col(2) + k.col(3)) * h / 6;

//   return ret;
// }

/*
// Qiita algorythm (1 step of Runge Kutta Fehlberg)
// k.col(0) = h * func(t, ret, ds);
// temp = ret + (k.col(0) * 0.25);
// k.col(1) = h * func(t + h * 0.25, temp, ds);
// temp = ret + ((k.col(0) * (0.09375)) + (k.col(1) * (0.28125)));
// k.col(2) = h * func(t + h * (0.375), temp, ds);
// temp = ret + (((0.87938) * k.col(0)) - ((3.27719) * k.col(1)) +
//               ((3.320892) * k.col(2)));
// k.col(3) = h * func(t + h * (0.92307), temp, ds);
// temp = ret + (((2.032407) * k.col(0)) - (8 * k.col(1)) +
//               ((7.173489) * k.col(2)) - ((0.20589) * k.col(3)));
// k.col(4) = h * func(t + h, temp, ds);
// temp = ret +
//        (-((0.29629) * k.col(0)) + (2 * k.col(1)) - ((1.38167) * k.col(2))
//        +
//         ((0.45330) * k.col(3)) - ((0.275) * k.col(4)));
// k.col(5) = h * func(t + h * 0.5, temp, ds);

// x_4th = ((0.11574 * k.col(0)) + (0.54892 * k.col(2)) +
//          (0.53533 * k.col(3)) - (0.2) * k.col(4));
// x_5th = ((0.1185 * k.col(0)) + (0.51898 * k.col(2)) + (0.50613 *
// k.col(3)) -
//          (0.18 * k.col(4)) + (0.03636 * k.col(5)));

// R = (x_5th - x_4th).norm() / h;
*/