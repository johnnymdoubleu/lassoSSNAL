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
SEXP transVecCpp(const Eigen::Map<Eigen::VectorXd> A){
  Eigen::setNbThreads(4);
  Eigen::VectorXd C = A.transpose();
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP crossCpp(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::setNbThreads(4);
  const int m(A.rows()), n(A.cols());
  Eigen::MatrixXd AtA(Eigen::MatrixXd(n, n).setZero().
                 selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()));
  
  return Rcpp::wrap(AtA);
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
SEXP eigenSparseMult(const Eigen::MappedSparseMatrix<double> A){
  
  Eigen::setNbThreads(4);
  Eigen::SparseMatrix<double> C = A.adjoint() * A;
  
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP eigenTransMapMatMult(const Eigen::Map<Eigen::MatrixXd> A,
                          int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A.transpose() * A;
  return Rcpp::wrap(C);
}