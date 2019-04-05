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
  double f0,f1,p0,p1,a,pran,pbase,ptop,prm,pout ;
  int selind = 0 ;

  NumericVector absc = df[0], f = df[1], pvec, pveccum ;

  for (int i = 0 ; i < lvec.size() ; i++) {
    pvec.push_back(lvec[i]/sum(lvec)) ;
    if (i > 0) {
      pveccum.push_back(pveccum[i - 1] + pvec[i]) ;
    } else pveccum.push_back(pvec[i]) ;
  }
  
  pveccum.push_front(0) ;
  
  pran = R::runif(0,1) ;
  
  for (int i = 0 ; i < pveccum.size() ; i++ ) {
    if (pveccum[i] <= pran) {
      selind = i ;
    }
  }
  
  pbase = pveccum[selind] ;
  ptop = pveccum[selind + 1] ;
  prm = pran - pbase ;
  
  if (absc[selind + 1] == zval && zind==1 && maxind > 0) {
    f0 = f[selind - 1] ;
    f1 = f[selind] ;
    p0 = absc[selind - 1] ;
    p1 = absc[selind] ;
    a = (f1 - f0)/(p1 - p0) ;
    pout = (log(exp(f1) + a*prm*sum(lvec)) + a*p1 - f1)/a ;
  } else if (absc[selind] == zval && zind==0 && maxind > 0) {
    f0 = f[selind + 1] ;
    f1 = f[selind + 2] ;
    p0 = absc[selind + 1] ;
    p1 = absc[selind + 2] ;
    a = (f1 - f0)/(p1 - p0) ;
    pout = (log(exp(f0 + a*(zval - p0)) + prm*a*sum(lvec)) + a*p0 - f0)/a ;
  } else {
    pout = absc[selind] + (absc[selind + 1] - absc[selind])*prm/(ptop - pbase) ;
  } 
  return pout ;
}

  