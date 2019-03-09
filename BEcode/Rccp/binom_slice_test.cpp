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
  double f0, f1, f2, fnew, p1, p2, pran, z ;
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
      p1 = std::max(pran - delta, 0.) ;
      p2 = std::min(pran + delta, 1.) ;
      
      f1 = R::dbinom(x, n, p1, false) ;
      f2 = R::dbinom(x, n, p2, false) ;
      
      while ( f1 > z && p1 > 0 ) {
        p1 = std::max(p1 - delta, 0.) ;
        f1 = R::dbinom(x, n, p1, false) ;
      }
      while ( f2 > z && p2 < 1) {
        p2 = std::min(p2 + delta, 1.) ;
        f2 = R::dbinom(x, n, p2, false) ;
      }
      fnew = 0;
      while ( fnew  < z ) {
        pran = R::runif(p1,p2) ;
        fnew = R::dbinom(x, n, pran, false) ;
      }
      pout[i] = pran ;
    }
    return pout ;
  }
}


/*** R
binomslicep(20,100,.3,.1,100)
*/
