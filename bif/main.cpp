#include "dynamical_system.hpp"
#include "newton.hpp"
#include <nlohmann/json.hpp>

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

  std::cout << "**************************************************"
            << std::endl;
  std::cout << "<initial condition>" << std::endl;
  std::cout << "ode solver : ";
  if (ds.use_classic_rk != false) {
    std::cout << "classic runge-kutta 4(5)" << std::endl;
  } else {
    std::cout << "runge-kutta-dormand-prince 5(6)" << std::endl;
  }
  std::cout << "increment parameter : " << ds.inc_param << std::endl;
  std::cout << "system dimention : " << ds.xdim << std::endl;
  std::cout << "period : " << ds.period << std::endl;
  std::cout << "p_sec def. : q(x) = x(" << ds.p_index << ") - " << ds.p_place
            << std::endl;
  std::cout << "params  : ";
  std::cout << ds.p.transpose() << std::endl;
  std::cout << "x0  : ";
  std::cout << ds.x0.transpose() << std::endl;
  std::cout << "tau : ";
  std::cout << ds.tau << std::endl;
  std::cout << "theta : ";
  std::cout << ds.theta << std::endl;
  std::cout << "**************************************************"
            << std::endl;

  newton(ds);

  // output json for latest state
  json["x0"] = ds.x0;
  json["params"] = ds.p;
  json["theta"] = ds.theta;
  std::ofstream json_out;
  json_out.open(ds.json_out_path, std::ios::out);
  json_out << json.dump(4);
  json_out.close();

  return 0;
}