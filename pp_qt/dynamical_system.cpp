#include "dynamical_system.hpp"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

Eigen::VectorXd dynamical_system::func(double t, const Eigen::VectorXd &x) {
  Eigen::VectorXd ret(xdim);

  // rossler
  ret(0) = -x(1) - x(2);
  ret(1) = x(0) + params(0) * x(1);
  ret(2) = params(1) * x(0) - params(2) * x(2) + x(0) * x(2);

  return ret;
}

double dynamical_system::q(const Eigen::VectorXd &x) {
  return x(p_index) - p_place;
}

dynamical_system::dynamical_system(const std::string &json_location) {
  std::ifstream ifs(json_location);
  if (ifs.fail()) {
    std::cerr << "File does NOT exist." << std::endl;
    std::exit(1);
  }

  nlohmann::json json;
  ifs >> json;

  xdim = json["x0"].size();

  p_index = json["p_index"];
  p_place = json["p_place"];
  direction = json["direction"];
  period = json["period"];

  tick = json["tick"];
  axis = json["axis"].get<std::vector<unsigned int>>();
  xrange = json["xrange"].get<std::vector<double>>();
  yrange = json["yrange"].get<std::vector<double>>();
  dparams = json["dparams"].get<std::vector<double>>();

  max_plot = json["max_plot"];
  max_poincare_plot = json["max_poincare_plot"];
  use_classic_rk = json["use_classic_rk"];
  rk_div = json["rk_div"];
  rkf_first_h = json["rkf_first_h"];
  rkf_h_max = json["rkf_h_max"];
  rkf_h_min = json["rkf_h_min"];
  rkf_tol = json["rkf_tol"];
  rkf_false_iter = json["rkf_false_iter"];
  poincare_eps = json["poincare_eps"];

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["x0"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;
  this->last_state = x0;

  std::vector<double> params_arr = json["params"];
  Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
  this->params = params;

  dqdx = Eigen::MatrixXd::Zero(1, xdim);
  dqdx(p_index) = 1.0;
}

void dynamical_system::integrate(double t0, const Eigen::VectorXd &x,
                                 double t_end) {
  Eigen::VectorXd state = x;
  Eigen::VectorXd next_state = x;
  Eigen::MatrixXd k = Eigen::MatrixXd::Zero(state.rows(), 6);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(state.rows());
  double t = t0;
  double h = rkf_first_h;
  double h_max = rkf_h_max;
  double h_min = rkf_h_min;
  double tol = rkf_tol;
  double R = 0.0;
  double delta = 0.0;
  unsigned int loop_counter = 0;
  bool flag = true;
  hit_section = false;

  while (flag) {
    k.col(0) = h * func(t0, state);
    temp = state + k.col(0) * 0.25;
    k.col(1) = h * func(t0 + h * 0.25, temp);
    temp = state + 0.09375 * k.col(0) + 0.28125 * k.col(1);
    k.col(2) = h * func(t0 + 0.375 * h, temp);
    temp = state + 0.8793809740555303 * k.col(0) -
           3.277196176604461 * k.col(1) + 3.3208921256258535 * k.col(2);
    k.col(3) = h * func(t0 + 0.9230769230769231 * h, temp);
    temp = state + 2.0324074074074074 * k.col(0) - 8.0 * k.col(1) +
           7.173489278752436 * k.col(2) - 0.20589668615984405 * k.col(3);
    k.col(4) = h * func(t0 + h, temp);
    temp = state - 0.2962962962962963 * k.col(0) + 2.0 * k.col(1) -
           1.3816764132553607 * k.col(2) + 0.4529727095516569 * k.col(3) -
           0.275 * k.col(3);
    k.col(5) = h * func(t0 + 0.5 * h, temp);

    R = (0.002777777777777778 * k.col(0) - 0.02994152046783626 * k.col(2) -
         0.029199893673577886 * k.col(3) + 0.02 * k.col(4) +
         0.03636363636363636 * k.col(5))
            .norm() /
        h;

    if (R <= tol) {
      loop_counter = 0;
      next_state += 0.11574074074074074 * k.col(0) +
                    0.5489278752436647 * k.col(2) +
                    0.5353313840155945 * k.col(3) - 0.2 * k.col(4);
      double qprod = q(next_state) * q(state);
      if (hit_section) {
        if (std::abs(q(next_state)) < poincare_eps) {
          state = next_state;
          x0 = state;
          t += h;
          QCPCsol.append(QCPCurveData(t, state(axis[0]), state(axis[1])));
          flag = false;
        }
        // Newton's method for h
        h -= (q(next_state) / (dqdx * func(h, state))(0, 0));
        next_state = state;
        continue;
      }
      if (qprod < 0 &&
          (direction < 0 ? next_state(p_index) - state(p_index)
                         : state(p_index) - next_state(p_index)) < 0 &&
          (next_state - state).norm() > 0.1) {
        hit_section = true;
        next_state = state;
        continue;
      }
      state = next_state;
      t += h;
      QCPCsol.append(QCPCurveData(t, state(axis[0]), state(axis[1])));
    }

    if (loop_counter++ > rkf_false_iter) {
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

    if (t >= t_end) {
      flag = false;
    } else if (t + h > t_end) {
      h = t_end - t;
      if (h < h_min) {
        flag = false;
      }
    }
  }
  last_state = state;
}

bool dynamical_system::is_hit_section() {
  if (hit_section) {
    hit_section = false;
    return true;
  } else {
    return false;
  }
}

void dynamical_system::integrate_rk45(double t0, const Eigen::VectorXd &x,
                                      double t_end) {
  Eigen::VectorXd state = x;
  Eigen::VectorXd next_state = x;
  Eigen::MatrixXd k(x.rows(), 4);
  Eigen::VectorXd temp(x.rows());
  double t = t0;
  double h = (t_end - t0) / rk_div;
  hit_section = false;

  for (int i = 0; i < rk_div; i++) {
    k.col(0) = func(t0, state);
    temp = state + h * 0.5 * k.col(0);
    k.col(1) = func(t0 + h * 0.5, temp);
    temp = state + h * 0.5 * k.col(1);
    k.col(2) = func(t0 + h * 0.5, temp);
    temp = state + h * k.col(2);
    k.col(3) = func(t0 + h, temp);

    next_state +=
        (h / 6.0) * (k.col(0) + 2.0 * k.col(1) + 2.0 * k.col(2) + k.col(3));
    double qprod = q(next_state) * q(state);
    if (hit_section) {
      if (std::abs(q(next_state)) < poincare_eps) {
        state = next_state;
        x0 = state;
        t += h;
        QCPCsol.append(QCPCurveData(t, state(axis[0]), state(axis[1])));
        break;
      }
      // Newton's method for h
      h -= (q(next_state) / (dqdx * func(h, state))(0, 0));
      next_state = state;
      continue;
    }
    if (qprod < 0 &&
        (direction < 0 ? next_state(p_index) - state(p_index)
                       : state(p_index) - next_state(p_index)) < 0 &&
        (next_state - state).norm() > 1.0e-02) {
      hit_section = true;
      next_state = state;
      continue;
    }
    state = next_state;
    t += h;
    QCPCsol.append(QCPCurveData(t, state(axis[0]), state(axis[1])));
  }
  last_state = state;
}