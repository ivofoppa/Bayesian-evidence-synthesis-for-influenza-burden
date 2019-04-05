#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double fpsample(DataFrame df, IntegerVector & parms, std::string dist, double zval, int maxind, int zind, NumericVector lvec) {
  double f0,f1,p0,p1,a ;

  NumericVector absc = df[0], f = df[1], pvec, pveccum ;

  pvec = lvec/sum(lvec) ;
  
  for 
  pveccum = c(0,cumsum(pvec))
  
  pran = runif(1)
  k = max(which(pveccum <= pran))
  pbase = pveccum[k]
  ptop = pveccum[k + 1]
  prm = pran - pbase
  
  if (absc[k + 1] == zval && zind==1 && maxind > 0) {
    f0 = f[k - 1]
    f1 = f[k]
    p0 = absc[k - 1]
    p1 = absc[k]
    a = (f1 - f0)/(p1 - p0)
    pout = (log(exp(f1) + a*prm*sum(lvec)) + a*p1 - f1)/a
  } else if (absc[k] == zval && zind==0 && maxind > 0) {
    f0 = f[k + 1]
    f1 = f[k + 2]
    p0 = absc[k + 1]
    p1 = absc[k + 2]
    a = (f1 - f0)/(p1 - p0)
    pout = (log(exp(f0 + a*(zval - p0)) + prm*a*sum(lvec)) + a*p0 - f0)/a
  } else {
    pout = absc[k] + (absc[k + 1] - absc[k])*prm/(ptop - pbase)
  } 
  pout
}

  