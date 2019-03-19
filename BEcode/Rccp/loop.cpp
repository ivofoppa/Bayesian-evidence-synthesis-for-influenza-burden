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

NumericVector samlingLoop(NumericVector absc,NumericVector lvec,NumericVector f,double zval,double maxind) {

  NumericVector abscaug  
  
}
  nsim <- 10000
    probls <- NULL 
    
    while (length(probls) < nsim) {
      psample <- fpsample(abscaug,lvec,f,zval,maxind)
      aran <- runif(1)
      lhull <- flowerhull(psample,abscaug,f)
      uhull <- fupperhull(psample,abscaug,f,zval,maxind)
      if (aran < lhull/uhull) {
        probls <- c(probls,psample)
      } else {
        aratio <- exp(fbin(psample))/uhull
        if (aran < aratio) {
          probls <- c(probls,psample)
        } else {
          abscaug <- unique(sort(round(c(abscaug,psample),digits = 8)))
          f <- sapply(abscaug,fbin)
          
          fmxind <- which(f==max(f))
          
          if (f[fmxind] > fbin(abscaug[fmxind] - 1e-10) & f[fmxind] > fbin(abscaug[fmxind] + 1e-10)) {
            maxind <- 0
          } else if (f[fmxind] < fbin(abscaug[fmxind] + 1e-5)) {
            maxind <- 1
          } else if (fbin(abscaug[fmxind] - 1e-5) > f[fmxind]) {
            maxind <- 2
          }
          
### calculating intersection x for middle section (containing max); using Rcpp function intsct2
          if (maxind==2) {
            zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
          } else {
            zabsc <- abscaug[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
          }
          
          zval <- intsct2(zabsc,x,n)
            
            abscaug <- sort(unique(c(abscaug,zval)))
            f <- sapply(abscaug,fbin)
            lvec <- fliksum(abscaug,f,zval,maxind)
        }
      }
    }
    
    mean(probls)
      hist(probls)
      

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
