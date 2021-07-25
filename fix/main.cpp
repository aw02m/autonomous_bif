#include "sys_common.hpp"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Put a path to the input json." << std::endl;
    std::exit(1);
  }

  std::ifstream ifs(argv[1]);
  if (ifs.fail()) {
    std::cerr << "File does NOT exist." << std::endl;
    std::exit(1);
  }

  nlohmann::json json;
  ifs >> json;
  dynamical_system ds(json);

  std::cout << "********************************" << std::endl;
  std::cout << "<Initial condition>" << std::endl;
  std::cout << "ODE solver : ";
  if (ds.use_classic_rk != false) {
    std::cout << "Classic Runge-Kutta 4(5)" << std::endl;
  } else {
    std::cout << "Runge-Kutta-Fehlberg 5(6)" << std::endl;
  }
  std::cout << "Increment parameter : " << ds.inc_param << std::endl;
  std::cout << "System dimention : " << ds.xdim << std::endl;
  std::cout << "period : " << ds.period << std::endl;
  std::cout << "p_sec def. : q(x) = x(" << ds.p_index << ") - " << ds.p_place
            << std::endl;
  std::cout << "params  : ";
  std::cout << ds.params.transpose() << std::endl;
  std::cout << "x0  : ";
  std::cout << ds.x0.transpose() << std::endl;
  std::cout << "u0  : ";
  std::cout << ds.u0.transpose() << std::endl;
  std::cout << "tau : ";
  std::cout << ds.tauk.transpose() << std::endl;
  std::cout << "********************************" << std::endl;

  newton(ds);

  return 0;
}