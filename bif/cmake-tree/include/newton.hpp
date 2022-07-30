#pragma once

#include "sys_common.hpp"
#include <chrono>

class dynamical_system;

void newton(dynamical_system &ds);
bool check_norm(Eigen::VectorXd F, dynamical_system &ds);