#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXcd expm_cplx(Eigen::MatrixXcd M) {
  return M.exp();
}