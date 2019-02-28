#include <Rcpp.h>

using namespace Rcpp;

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
    f = R::choose(n,x)*pow(p,x)*pow((1-p),(n-x)) ;
    return(f) ;
  }
}
