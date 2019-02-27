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
int choosecpp(int n, int k) { 
  // Base Cases 
  if (k==0 || k==n) 
    return 1; 
  if ( k > (n - k))
    k = n - k ;
  // Recur 
  return  choosecpp(n-1, k-1) + choosecpp(n-1, k); 
} 
// [[Rcpp::export]]
double binompdf(int x, int n, double p) 
{
  double f ;
  if ( (p < 0) || (p > 1) || (x > n) ) 
  {
    return(0) ; 
  }
  else 
  {
    f = choosecpp(n,x)*pow(p,x)*pow((1-p),(n-x)) ;
    return(f) ;
  }
}
