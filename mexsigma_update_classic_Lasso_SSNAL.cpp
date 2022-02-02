#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

#if !defined(max)
#define max(a,b) (a>b?a:b)
#endif

#if !defined(min)
#define min(a,b) (a<b?a:b)
#endif

// [[Rcpp::export]]
NumericVector mexsigma_update(double sigma, double sigmamax, double sigmamin, double prim_win, double dual_win,
                               int iter, int inner_breakyes) {
  /****** main body ****************/ 
  int sigma_update_iter;
  double sigmascale = 5.0;
  
  if ( iter < 10)
    sigma_update_iter = 2;
  else if(iter < 20)
    sigma_update_iter = 3;
  else if(iter < 200)
    sigma_update_iter = 3;
  else if(iter < 500)
    sigma_update_iter = 10;
  else
    sigma_update_iter = 20;
  
  if ((iter%sigma_update_iter == 0) && (inner_breakyes < 0)) {
    if (prim_win > max(1, 1.2*dual_win)){
      prim_win = 0.0;
      sigma = min(sigmamax,sigma*sigmascale);
    }else if (dual_win > max(1, 1.2*prim_win)){
      dual_win = 0.0;
      sigma = max(sigmamin,sigma/sigmascale);
    } 
  }
  
  if (inner_breakyes >= 0 && iter >= 10) {
    sigma = max(sigmamin,2*sigma/sigmascale);
  }
  
  NumericVector out(3);
  
  out[0]=sigma;
  out[1]=prim_win;
  out[2]=dual_win;
  
  return out;
}

// [[Rcpp::export]]
NumericVector mexscale(double sigma, double normx, double normuxi,
                       double bscale, double cscale) {
  double bscale2, cscale2, cst = 1.0;
  double sbc, sboc;
  
  /****** main body ****************/ 
  if (normx < pow(10,-7)) { normx = 1; normuxi = 1; }
  bscale2 = normx*cst;
  cscale2 = normuxi*cst;
  sbc  = sqrt(bscale2*cscale2);       
  sboc = sqrt(bscale2/cscale2);
  sigma = sigma*(cscale2/bscale2);
  cscale = cscale*cscale2;
  bscale = bscale*bscale2;
  
  NumericVector out(7);
  
  out[0] = sigma;
  out[1] = bscale2;
  out[2] = cscale2;
  out[3] = sbc;
  out[4] = sboc;
  out[5] = bscale;
  out[6] = cscale;
  
  return out;
}