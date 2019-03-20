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
// [[Rcpp::export]]
void rcpp_rprintf(IntegerVector v){
  // printing values of all the elements of Rcpp vector  
  for(int i=0; i<v.length(); ++i){
    Rprintf("the value of v[%i] : %i \n", i, v[i]);
  }
}