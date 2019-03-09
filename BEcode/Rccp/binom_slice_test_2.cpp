#include <Rcpp.h>
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
NumericVector binomslicep(int x, int n, double p0, double delta, int num ) 
{
  double f0, logf1, logf2, logfnew, p1, p2, pran, z, logz ;
  NumericVector pout(num) ;
  
  if ( (p0 < 0) || (p0 > 1) || (x > n)) 
  {
    return(0) ; 
  }
  else 
  {
    pran = p0 ;
    for ( int i=0 ; i <= num ; i++ ) {
      f0 = R::dbinom(x, n, pran, false) ;
      z = R::runif(0,f0) ;
      logz = std::log(z) ;
      p1 = std::max(pran - delta, 0.) ;
      p2 = std::min(pran + delta, 1.) ;
      
      logf1 = R::dbinom(x, n, p1, true) ;
      logf2 = R::dbinom(x, n, p2, true) ;
      
      while ( logf1 > logz && p1 > 0 ) {
        p1 = std::max(p1 - delta, 0.) ;
        logf1 = R::dbinom(x, n, p1, true) ;
      }
      while ( logf2 > logz && p2 < 1) {
        p2 = std::min(p2 + delta, 1.) ;
        logf2 = R::dbinom(x, n, p2, true) ;
      }
      logfnew = -9999;
      while ( logfnew  < logz ) {
        pran = R::runif(p1,p2) ;
        logfnew = R::dbinom(x, n, pran, true) ;
      }
      pout[i] = pran ;
    }
    return pout ;
  }
}


/*** R
binomslicep(20,100,.3,.1,100)
*/
