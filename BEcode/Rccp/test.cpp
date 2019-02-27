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
int choosecpp(int n, int k)
{
  { 
    int res = 1; 
    
    // Since C(n, k) = C(n, n-k) 
    if ( k > n - k ) 
      k = n - k; 
    
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1] 
    for (int i = 0; i < k; i++) 
    { 
      res *= (n - i); 
      res /= (i + 1); 
    } 
    
    return res; 
  } 
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
