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
int abscPrep(NumericVector absc, NumericVector f, int x, int n ){
  int maxind = -99 ;
  int fmxind = which_max(f) ;
  double zabsc ;
  double arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
  
  if (f[fmxind] > R::dbinom(x,n,arg1,true) && f[fmxind] > R::dbinom(x,n,arg2,true)) {
    maxind = 0 ;
  } else if (f[fmxind] < R::dbinom(x,n,arg2,true)) {
    maxind = 1 ;
  } else if (R::dbinom(x,n,arg1,true) > f[fmxind]) {
    maxind = 2 ;
    }

  int minpexp = -10 ;
//* calculating intersection x for middle section (containing max); using Rcpp function intsct2
  if  (maxind==2 ) {
    if ( fmxind > 1 && fmxind < (absc.size() - 1)) {
      zabsc = absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
    } else if (fmxind <= 2) {
      while (fmxind <= 2) {
        eval(parse(text=paste0('minp <- 1e',minpexp)))
        absc <- sort(unique(c(seq(minp,absc[2],length.out = 4),absc)))
        f <- sapply(absc,fbin)
        fmxind <- which(f==max(f))
        minpexp <- minpexp*2
      }
      zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
    } else if (fmxind==length(absc)) {
      while (fmxind==length(absc)) {
        absc <- sort(unique(c(seq(tail(absc,2)[1],tail(absc,1),length.out = 3),absc)))
        f <- sapply(absc,fbin)
        fmxind <- which(f==max(f))
      }
      zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
    }
  } else {
    if (fmxind > 1 & fmxind < (length(absc) - 1)) {
      zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
    } else if (fmxind == 1) {
      while (fmxind==1) {
        absc <- sort(unique(c(seq(absc[1],absc[2],length.out = 3),absc)))
        f <- sapply(absc,fbin)
        fmxind <- which(f==max(f))
      }
      zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
    } else if (fmxind > (length(absc) - 2)) {
      while (fmxind > (length(absc) - 2)) {
        absc <- sort(unique(c(seq(tail(absc,2)[1],tail(absc,1),length.out = 4),absc)))
        f <- sapply(absc,fbin)
        fmxind <- which(f==max(f))
      }
      zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
    }
  }
  
  
  
  
      return maxind ;
}


/*** R
absc <- c(1e-40,.05,.1,.15,.21,.25,.3,.4,.5,.6,.7,.8,.9,1 - 1e-40)

f <- sapply(absc,function(p) dbinom(30,100,p,TRUE))

abscPrep(absc,f,30,100)
*/

