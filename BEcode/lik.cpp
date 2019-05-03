#include <Rcpp.h>
#include <cmath>
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
double lik(double pr, List & parms, std::string liktype) {
  double L = 999 ;
  if (liktype=="ptest") {
    NumericVector ptest = parms["ptest"], ntot = parms["ntot"], sens = parm["sens"], pvec ;
    IntegerVector  ttest = parms["ttest"];
    int ind = parms["ind"] ;
    
    ptest[ind-1] = int (pr) ;
    
    for (int k=0 ; k < sens.size() ; k++) {
      pvec.push_back(sens[k] * ptest[k]) ;
    }
    
    double pt = sum(pvec) ;
    int x1 = ttest[ind-1], N1 = ntot[ind-1] ;
    pvec = NumericVector::create() ;
    double pmultinom = 1;
    for (int k = 0; k < ptest.size() ; k++) {
      pmultinom *= pow(ptest[k],ttest[k]) ;
    }
    L = R::dbinom(x1,N,pr,0) * pmultinom ;

      } else if (liktype=="binomN") {
    int x = parms["x"];
        
    double p = parms["p"] ;
    
    L = R::dbinom(x,pr,p,0) ;
  }
  return L ;
}
/***R
parms <- list(x=20,N=100,p=.25)

lik(.3,parms,"binomp")
dbinom(20,100,.3)

lik(100,parms,"binomN")
dbinom(20,100,.25)
*/