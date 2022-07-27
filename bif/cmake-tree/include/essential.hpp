#include <eigen3/Eigen/Dense>

double det_derivative(const Eigen::MatrixXd &A, const Eigen::MatrixXd &dA,
                      unsigned int dim) {
  Eigen::MatrixXd temp(dim, dim);
  double ret = 0;

  for (int i = 0; i < dim; i++) {
    temp = A;
    temp.col(i) = dA.col(i);
    ret += temp.determinant();
  }

  return ret;
}

Eigen::MatrixXd bialt_prod(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                           unsigned int xdim, unsigned int bialt_dim) {
  Eigen::MatrixXd ret(bialt_dim, bialt_dim);

  int row = 0;
  int col = 0;
  Eigen::Matrix2d temp;
  double aug_det;
  for (int p = 1; p < xdim; p++) {
    for (int q = 0; q < p; q++) {
      for (int r = 1; r < xdim; r++) {
        for (int s = 0; s < r; s++) {
          temp(0, 0) = A(p, r);
          temp(0, 1) = A(p, s);
          temp(1, 0) = B(q, r);
          temp(1, 1) = B(q, s);
          aug_det = temp.determinant();
          temp(0, 0) = B(p, r);
          temp(0, 1) = B(p, s);
          temp(1, 0) = A(q, r);
          temp(1, 1) = A(q, s);
          aug_det += temp.determinant();
          aug_det /= 2;
          ret(row, col++) = aug_det;
        }
      }
      col = 0;
      row++;
    }
  }

  return ret;
}

Eigen::MatrixXd biproduct(const Eigen::MatrixXd &A, unsigned int xdim,
                          unsigned int bialt_dim) {
  Eigen::MatrixXd ret(bialt_dim, bialt_dim);

  int row = 0;
  int col = 0;
  for (int p = 1; p < xdim; p++) {
    for (int q = 0; q < p; q++) {
      for (int r = 1; r < xdim; r++) {
        for (int s = 0; s < r; s++) {
          if (r == q) {
            ret(row, col++) = -A(p, s);
          } else if (r != p && s == q) {
            ret(row, col++) = A(p, r);
          } else if (r == p && s == q) {
            ret(row, col++) = A(p, p) + A(q, q);
          } else if (r == p && s != q) {
            ret(row, col++) = A(q, s);
          } else if (s == p) {
            ret(row, col++) = -A(q, r);
          } else {
            ret(row, col++) = 0;
          }
        }
      }
      col = 0;
      row++;
    }
  }

  return ret;
}

Eigen::MatrixXd biproduct_derivative(const Eigen::MatrixXd &dA,
                                     unsigned int xdim,
                                     unsigned int bialt_dim) {
  Eigen::MatrixXd ret(bialt_dim, bialt_dim);

  int row = 0;
  int col = 0;
  for (int p = 1; p < xdim; p++) {
    for (int q = 0; q < p; q++) {
      for (int r = 1; r < xdim; r++) {
        for (int s = 0; s < r; s++) {
          if (r == q) {
            ret(row, col++) = -dA(p, s);
          } else if (r != p && s == q) {
            ret(row, col++) = dA(p, r);
          } else if (r == p && s == q) {
            ret(row, col++) = dA(p, p) + dA(q, q);
          } else if (r == p && s != q) {
            ret(row, col++) = dA(q, s);
          } else if (s == p) {
            ret(row, col++) = -dA(q, r);
          } else {
            ret(row, col++) = 0;
          }
        }
      }
      col = 0;
      row++;
    }
  }

  return ret;
}

Eigen::MatrixXd bialt_prod_square(const Eigen::MatrixXd &A, unsigned int xdim,
                                  unsigned int bialt_dim) {
  Eigen::MatrixXd ret(bialt_dim, bialt_dim);

  int row = 0;
  int col = 0;
  Eigen::Matrix2d temp;
  for (int p = 1; p < xdim; p++) {
    for (int q = 0; q < p; q++) {
      for (int r = 1; r < xdim; r++) {
        for (int s = 0; s < r; s++) {
          temp(0, 0) = A(p, r);
          temp(0, 1) = A(p, s);
          temp(1, 0) = A(q, r);
          temp(1, 1) = A(q, s);
          ret(row, col++) = temp.determinant();
        }
      }
      col = 0;
      row++;
    }
  }

  return ret;
}

Eigen::MatrixXd bialt_prod_square_derivative(const Eigen::MatrixXd &A,
                                             const Eigen::MatrixXd &dA,
                                             unsigned int xdim,
                                             unsigned int bialt_dim) {
  Eigen::MatrixXd ret(bialt_dim, bialt_dim);

  int row = 0;
  int col = 0;
  Eigen::Matrix2d temp;
  Eigen::Matrix2d dtemp;
  for (int p = 1; p < xdim; p++) {
    for (int q = 0; q < p; q++) {
      for (int r = 1; r < xdim; r++) {
        for (int s = 0; s < r; s++) {
          temp(0, 0) = A(p, r);
          temp(0, 1) = A(p, s);
          temp(1, 0) = A(q, r);
          temp(1, 1) = A(q, s);
          dtemp(0, 0) = dA(p, r);
          dtemp(0, 1) = dA(p, s);
          dtemp(1, 0) = dA(q, r);
          dtemp(1, 1) = dA(q, s);
          ret(row, col++) = det_derivative(temp, dtemp, 2);
        }
      }
      col = 0;
      row++;
    }
  }

  return ret;
}

unsigned int nearest_idx(const Eigen::VectorXcd &array, double target_val) {
  unsigned int target_index = 0;

  double delta = std::abs(array(0).real()) - target_val;
  for (int i = 1; i < array.rows(); i++) {
    if (std::abs(array(i).real()) < delta) {
      target_index = i;
      delta = std::abs(array(i).real()) - target_val;
    }
  }

  return target_index;
}