#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector binomslicep(int x, int n, double p0, double delta, int num ) 
{
  double logf0, logf1, logf2, logfnew, yran, p1, p2, pran ;
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
      yran = R::rexp(1) ;
      p1 = std::max(p0 - delta, 0.) ;
      p2 = std::min(p0 + delta, 1.) ;
      
      f1 = R::dbinom(x, n, p1, true) ;
      f2 = R::dbinom(x, n, p2, true) ;
      
      while ( f1 > yran && p1 > 0 ) {
        p1 = std::max(p1 - delta, 0.) ;
        f1 = R::dbinom(x, n, p1, true) ;
      }
      while ( f2 > yran && p2 < 1) {
        p2 = std::min(p2 + delta, 1.) ;
        f2 = R::dbinom(x, n, p2, true) ;
      }
      fnew = 0.;
      while ( fnew  < yran ) {
        pran = R::runif(p1,p2) ;
        fnew = R::dbinom(x, n, pran, true) ;
      }
      pout[i] = pran ;
    }
    return pout ;
  }
}
// [[Rcpp::export]]
IntegerVector binomslicen(int x, int n0, double p, int delta, int num ) 
{
  double f0, f1, f2, fnew, yran, numran ;
  int n1, n2, nran ;
  IntegerVector nout(num) ;
  
  if ( (p < 0) || (p > 1) || (n0 < 1) || (x < 0) || (n0 < x)) 
  {
    return(0) ; 
  }
  else 
  {
    nran = n0 ;
    for (int i=0 ; i <= num ; i++ ) {
      f0 = R::dbinom(x, nran, p, log) ;
      yran = R::runif(0,f0) ;
      n1 = std::max(nran - delta, 0) ;
      n2 = nran + delta ;
      
      f1 = R::dbinom(x, n1, p, false) ;
      f2 = R::dbinom(x, n2, p, false) ;
      
      while ( f1 > yran && (n1 - delta) >= x ) {
        n1 = std::max(n1 - delta, x) ;
        f1 = R::dbinom(x, n1, p, false) ;
      }
      n1 = std::max(n1, x) ;
      while ( f2 > yran ) {
        n2 = n2 + delta ;
        f2 = R::dbinom(x, n2, p, false) ;
      }
      fnew = 0.;
      while ( fnew < yran ) {
        numran = R::runif(0,1) ;
        numran = numran * (n2 - n1) + n1 ;
        nran = round(numran) ;
        fnew = R::dbinom(x, nran, p, false) ;
      }
      nout[i] = nran ;
    }
    return(nout) ;
  }
}
