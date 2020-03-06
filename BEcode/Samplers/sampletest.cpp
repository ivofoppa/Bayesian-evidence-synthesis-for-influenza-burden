#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector samples(NumericVector & absc, IntegerVector & parms, std::string dist, int nsims) {
  int maxind = -99, zind = -99 ;
  double zval ;
  NumericVector zabsc, lvec, parmsmplvec ;
  DataFrame df ;
  
  Function abscPrep("abscPrep") ;
  Function maxindC("maxindC") ;
  Function zabscPrep("zabscPrep") ;
  Function intsct("intsct") ;
  Function fabscaug("fabscaug") ;
  Function fzind("fzind") ;
  Function fliksum("fliksum") ;
  Function fpsample("fpsample") ;
  
  
  while (parmsmplvec.size() < nsims ) {
    df = abscPrep(absc,parms,dist) ;
    
    maxind = maxindC(df,parms,dist) ;
    
    zabsc = zabscPrep(df,parms,dist,maxind) ;
    
    zval = intsct(zabsc,parms, dist) ;
    
    df = fabscaug(df,parms,dist,zval) ;
    
    zind = fzind(df,parms,dist,zval) ;
    lvec = fliksum(df, parms, dist, zval, maxind,zind) ;
    
    parmsmplvec.push_back(fpsample(df, parms, dist, zval, maxind,zind,lvec)) ;
  }
  return parmsmplvec ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
parms <- c(10,100)
dist <- "Poi"
dist <- "binom"
absc <- c(.05,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9)

*/
