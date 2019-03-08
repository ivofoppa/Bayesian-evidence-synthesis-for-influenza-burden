#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector binomslicep(int x, int n, double p0, double delta, int num ) 
{
  double logf0, logf1, logf2, logfnew, p1, p2, pran, z, e ;
  NumericVector pout(num) ;
  
  if ( (p0 < 0) || (p0 > 1) || (x > n)) 
  {
    return(0) ; 
  }
  else 
  {
    pran = p0 ;
    for ( int i=0 ; i <= num ; i++ ) {
      logf0 = R::dbinom(x, n, pran, true) ;
      e = R::rexp(1) ;
      z = logf0 + e ;
      p1 = std::max(p0 - delta, 0.) ;
      p2 = std::min(p0 + delta, 1.) ;
      
      logf1 = R::dbinom(x, n, p1, true) ;
      logf2 = R::dbinom(x, n, p2, true) ;
      
      while ( logf1 > z && p1 > 0 ) {
        p1 = std::max(p1 - delta, 0.) ;
        logf1 = R::dbinom(x, n, p1, true) ;
      }
      while ( logf2 > z && p2 < 1) {
        p2 = std::min(p2 + delta, 1.) ;
        logf2 = R::dbinom(x, n, p2, true) ;
      }
      logfnew = -99999;
      while ( logfnew  < z ) {
        pran = R::runif(p1,p2) ;
        logfnew = R::dbinom(x, n, pran, true) ;
      }
      pout[i] = pran ;
    }
    return pout ;
  }
}
