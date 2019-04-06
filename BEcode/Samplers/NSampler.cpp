#include <Rcpp.h>
#include <limits>
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
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
int binomnSampler(NumericVector parms2){
  int ntest0 = 0, ntest = 0, outn = 0;
  double x = parms2[0], p = parms2[1], ptest0 =0, ptest = 0 ;
  IntegerVector nvec ;
  NumericVector pvec ;
  
  ptest0 = R::dbinom(x,ntest0,p,true) ;
  ntest0 = round(x/p) ;
  ntest = ntest0 + 1 ;
  ptest = R::dbinom(x,ntest,p,true) ;
  
  pvec = NumericVector::create(exp(ptest0)) ;
  nvec = IntegerVector::create(ntest0) ;
  
  while (ptest > -50) {
    nvec.push_back(ntest) ;
    pvec.push_back(exp(ptest)) ;
    ntest++ ;
    ptest = R::dbinom(x,ntest,p,true) ;
  }
  
  ntest = ntest0 - 1 ;
  ptest= R::dbinom(x,ntest,p,true) ;
  
  while (ptest > -50 && ntest >= x) {
    nvec.push_front(ntest);
    pvec.push_front(exp(ptest)) ;
    ntest-- ;
    ptest = R::dbinom(x,ntest,p,true) ;
  }
  
  int nsz = nvec.size();
  double total = 0;
  for(int i = 0 ; i < nsz ; ++i) {
    total += pvec[i];
  }
  
  for(int i = 0 ; i < nsz ; ++i) {
    pvec[i] /= total ;
  }
  
  double uran =  R::runif(0,1) ;
  
  double psum = 0 ;
  int selint = 1; 
  int i = 0 ;
  while (selint == 1) {
    psum += pvec[i] ;
    if (psum >=  uran) {
      selint = 0 ;
    }
    i++ ;
  }
  
  outn = nvec[i - 2] ;
  return outn ;
}
/***R
parms <- c(10,.2)
binomnSampler(parms)
*/
