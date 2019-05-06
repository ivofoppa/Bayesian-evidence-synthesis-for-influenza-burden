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
double lik(double pr,List & parms,List & inits,std::string liktype) {
  double L = 999 ;
  if (liktype=="ptest1") {
    NumericVector ptest = parms["ptest"],sens = parms["sens"],pvec ;
    IntegerVector  nttype = parms["nttype"];
    int ind = parms["ind"],FSNfluhosp = parms["FSNfluhosp"],
                    fluhosp = parms["fluhosp"] ;
    double ptestnorm ;
    
    ptest[ind-1] = pr ;
    ptestnorm = sum(ptest) ;
    
    for (int k=0 ; k < sens.size() ; k++) {
      pvec.push_back(sens[k] * ptest[k]) ;
    }
    
    double pt = sum(pvec) ;
    pvec = NumericVector::create() ;
    double pmultinom = 0;
    for (int k = 0; k < ptest.size() ; k++) {
      pmultinom += nttype[k] * log(ptest[k]/ptestnorm) ;
    }
    L = R::dbinom(FSNfluhosp,fluhosp,pt,1) + pmultinom ;

      } else if (liktype=="binomN") {
    int x = parms["x"];
        
    double p = parms["p"] ;
    
    L = R::dbinom(x,pr,p,1) ;
  }
  return L ;
}
/***R
nttype <- nttype[1,]
ptest <- nttype/sum(nttype)

sens <- c(0.86,)
parms <- list(x=20,N=100,p=.25)

lik(100,parms,"binomN")
dbinom(20,100,.25,log = T)
*/