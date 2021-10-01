#ifndef SYS_COMMON_HPP_
#define SYS_COMMON_HPP_

#define EIGEN_NO_DEBUG

#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>

#define debug(var)                                                             \
  do {                                                                         \
    std::cout << #var << " : " << std::endl;                                   \
    view(var);                                                                 \
  } while (0)
#define debug_exit(var)                                                        \
  do {                                                                         \
    std::cout << #var << " : " << std::endl;                                   \
    view(var);                                                                 \
    exit(0);                                                                   \
  } while (0)
template <typename T> void view(T e) { std::cout << e << std::endl; }
template <typename T> void view(const std::vector<T> &v) {
  for (const auto &e : v) {
    std::cout << e << " ";
  }
  std::cout << std::endl;
}
template <typename T> void view(const std::vector<std::vector<T>> &vv) {
  for (const auto &v : vv) {
    view(v);
  }
}

#endif