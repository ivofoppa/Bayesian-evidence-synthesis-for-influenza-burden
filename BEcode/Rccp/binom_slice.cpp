#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double binompmf(int x, int n, double p) 
{
  double f ;
  if ( (p < 0) || (p > 1) || (x > n) ) 
  {
    return(0) ; 
  }
  else 
  {
    f = R::choose(n,x)*pow(p,x)*pow((1-p),(n-x)) ;
    return(f) ;
  }
}

// [[Rcpp::export]]
NumericVector binomslicep(int x, int n, double p0, double delta, int num ) 
{
  double f0, f1, f2, fnew, yran, p1, p2, pran ;
  NumericVector pout(num) ;
  
  if ( (p0 < 0) || (p0 > 1) || (x > n)) 
  {
    return(0) ; 
  }
  else 
  {
    pran = p0 ;
    for ( int i=0 ; i <= num ; i++ ) {
      f0 = binompmf(x, n, pran) ;
      yran = R::runif(0,f0) ;
      p1 = std::max(p0 - delta, 0.) ;
      p2 = std::min(p0 + delta, 1.) ;
      
      f1 = binompmf(x, n, p1) ;
      f2 = binompmf(x, n, p2) ;
      
      while ( f1 > yran && p1 > 0 ) {
        p1 = std::max(p1 - delta, 0.) ;
        f1 = binompmf(x, n, p1) ;
      }
      while ( f2 > yran && p2 < 1) {
        p2 = std::min(p2 + delta, 1.) ;
        f2 = binompmf(x, n, p2) ;
      }
      fnew = 0.;
      while ( fnew  < yran ) {
        pran = R::runif(p1,p2) ;
        fnew = binompmf(x, n, pran) ;
      }
      pout[i] = pran ;
    }
    return(pout) ;
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
      f0 = binompmf(x, nran, p) ;
      yran = R::runif(0,f0) ;
      n1 = std::max(nran - delta, 0) ;
      n2 = nran + delta ;
      
      f1 = binompmf(x, n1, p) ;
      f2 = binompmf(x, n2, p) ;
      
      while ( f1 > yran && (n1 - delta) >= x ) {
        n1 = std::max(n1 - delta, x) ;
        f1 = binompmf(x, n1, p) ;
      }
      n1 = std::max(n1, x) ;
      while ( f2 > yran ) {
        n2 = n2 + delta ;
        f2 = binompmf(x, n2, p) ;
      }
      fnew = 0.;
      while ( fnew < yran ) {
        numran = R::runif(0,1) ;
        numran = numran * (n2 - n1) + n1 ;
        nran = round(numran) ;
        fnew = binompmf(x, nran, p) ;
      }
      nout[i] = nran ;
    }
    return(nout) ;
  }
}
