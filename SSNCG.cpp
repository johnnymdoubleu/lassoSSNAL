// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenVecMatMult(Eigen::Map<Eigen::MatrixXd> A,
                     Eigen::Map<Eigen::VectorXd> b,
                     int n_cores){

}