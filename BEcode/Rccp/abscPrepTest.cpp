#include <Rcpp.h>
#include <limits>
#include <vector>
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
NumericVector abscPrep(NumericVector & absc, NumericVector & f, int x, int n ){
  int maxind = -99, fmxind = which_max(f), minpexp = -5 ;
  NumericVector zabsc(4) ;
  double arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
  double minp = 0, pintv = 0, abscinc = 0 ;
  
  if (f[fmxind] > R::dbinom(x,n,arg1,true) && f[fmxind] > R::dbinom(x,n,arg2,true)) {
    maxind = 0 ;
  } else if (f[fmxind] < R::dbinom(x,n,arg2,true)) {
    maxind = 1 ;
  } else if (R::dbinom(x,n,arg1,true) > f[fmxind]) {
    maxind = 2 ;
  }
  
  //* calculating intersection x for middle section (containing max); using Rcpp function intsct2
  if  (maxind==2 ) {
    if (fmxind < 2) {
      
      while (fmxind < 2) {
        minp = pow(10,minpexp) ;
        pintv = absc[1] - minp ;
        absc.erase(0) ;
        abscinc = pintv/2 ;
        absc.insert(absc.begin(), minp) ;
        absc.insert(absc.begin(), minp + abscinc) ;
        
        NumericVector fnew = (R::dbinom(x,n,absc[0],true)) ;
        for (int i=1 ; i < absc.size() ; i++) {
          fnew.push_back(R::dbinom(x,n,absc[i],true)) ;
        }
        fmxind = which_max(fnew) ;
        minpexp = minpexp*2 ;
      }
    } else if (fmxind==(absc.size() - 1)) {
      while (fmxind==(absc.size() - 1)) {
        pintv = 1 - absc[absc.size() - 1] ;
        abscinc = pintv/2 ;
        absc.push_back(absc[absc.size() - 1] + abscinc) ;
        absc.push_back(absc[absc.size() - 1] + abscinc) ;
        
        NumericVector fnew = (R::dbinom(x,n,absc[0],true)) ;
        for (int i=1; i < absc.size() ; i++) {
          fnew.push_back(R::dbinom(x,n,absc[i],true)) ;
        }
        fmxind = which_max(fnew) ;
      }
    }
    for (int i=0 ; i < 3 ; i++) {
      zabsc[0] = absc[fxmind - 2 + i] ;
    }
  } else {
    
    if (fmxind == 0) {
      while (fmxind==0) {
        minp = pow(10,minpexp) ;
        pintv = absc[1] - minp ;
        absc.erase(0) ;
        abscinc = pintv/2 ;
        absc.insert(absc.begin(), minp) ;
        absc.insert(absc.begin(), minp + abscinc) ;
        
        NumericVector fnew = (R::dbinom(x,n,absc[0],true)) ;
        for (int i=1 ; i < absc.size() ; i++) {
          fnew.push_back(R::dbinom(x,n,absc[i],true)) ;
        }
        fmxind = which_max(fnew) ;
        minpexp = minpexp*2 ;
      }
    } else if (fmxind > (absc.size() - 3)) {
      while (fmxind > (absc.size() - 3)) {
        absc <- sort(unique(c(seq(tail(absc,2)[1],tail(absc,1),length.out = 4),absc)))
        f <- sapply(absc,fbin)
        fmxind <- which(f==max(f))
      }
      zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
    }
  }
  
  zabsc = (absc[fmxind - 2],absc[fmxind - 1],absc[fmxind],absc[fmxind + 1]) ;
} 

/*** R
absc <- c(0.2,.05,.1,.15,.21,.25,.3,.4,.5,.6,.7,.8,.9,1 - 1e-40)

x <- 1 ; n <- 100
f <- sapply(absc,function(p) dbinom(x,n,p,TRUE))

abscPrep(absc,f,x,n)
absc
  */

