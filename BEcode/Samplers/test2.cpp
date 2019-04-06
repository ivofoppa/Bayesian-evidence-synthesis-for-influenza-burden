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
double fbin(double p, NumericVector parms){
  double x = parms[0], n = parms[1] ;
  double f = R::dbinom(x,n,p,false) ;
  return f ;
} 

// [[Rcpp::export]]
double fpoi(double r, NumericVector parms){
  double x = parms[0] ;
  double f = R::dpois(x,r,false) ;
  return f ;
} 

// [[Rcpp::export]]
SEXP f(double pr, NumericVector parms, std::string dist){
  Function fn1("fbin") ;
  Function fn2("fpoi") ;
  SEXP fv = 0 ;
  if (dist=="binom") {
    fv = fn1(pr,parms) ;
  } else if (dist=="Poi") {
    fv = fn2(pr,parms) ;
  }
  return fv ;
} 

/***R
parms <- c(10,100)
  #dist <- "Poi"
dist <- "binom"
 p <- .3
  
 f(.3,parms,dist)
  */
