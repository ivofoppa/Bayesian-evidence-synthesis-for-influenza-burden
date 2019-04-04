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
NumericVector zabsc2(NumericVector & absc, int x, int n ){
  Function f("zabsc") ;
  NumericVector zabsc = f(absc,x,n) ;
  return zabsc ;
} 

/***R
parms <- c(10,100)
  dist <- "Poi"
dist <- "binom"
absc <- c(.05,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9)
  
  df <- abscPrep(absc,parms,dist)
  
  maxind <- maxindC(df,parms,dist)
  
  zabsc <- zabscPrep(df,parms,dist,maxind)
  
  zval <- intsct(zabsc,parms, dist)
  
  f0 <- sapply(df[,1],function(p) dbinom(10,100,p,log = T))
  
  fabscaug(df,parms,dist,zval)
  
  flowerhull(.2,df,parms)
  
  f[4] > dbinom(1,100,.0125,log = T)
  */
