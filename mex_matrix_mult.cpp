// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenVecMatMult(Eigen::Map<Eigen::MatrixXd> A,
                     Eigen::Map<Eigen::VectorXd> b,
                     int n_cores){
  Eigen::setNbThreads(n_cores);
  Eigen::VectorXd C = A * b;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP transCpp(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::setNbThreads(4);
  Eigen::MatrixXd C = A.transpose();
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A,
                     const Eigen::Map<Eigen::MatrixXd> B,
                     int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP eigenTransMapMatMult(const Eigen::Map<Eigen::MatrixXd> A,
                          int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A.transpose() * A;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMult(const Eigen::Map<Eigen::MatrixXd> A){
  
  Eigen::setNbThreads(4);
  Eigen::MatrixXd C = A.transpose() * A;
  return Rcpp::wrap(C);
}