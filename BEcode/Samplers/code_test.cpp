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
  double fbin(double p, NumericVector parms){
    double x = parms[0], n = parms[1] ;
    double f = R::dbinom(x,n,p,false) ;
    return f ;
  } 
  
  double fpoi(double r, NumericVector parms){
    double x = parms[0] ;
    double f = R::dpois(x,r,false) ;
    return f ;
  } 
  
  double fdist(double pr, NumericVector parms, std::string dist){
    double fv = 0 ;
    if (dist=="binom") {
      fv = fbin(pr,parms) ;
    } else if (dist=="Poi") {
      fv = fpoi(pr,parms) ;
    }
    return fv ;
  } 
  // [[Rcpp::export]]
DataFrame abscPrep(NumericVector & absc, NumericVector & parms, std::string & dist ){
  int maxind = 0, fmxind = 0 ;
  NumericVector f ;
  double arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
  double minp = 0, pintv = 0, abscinc = 0, newabsc = 0 ;
  
  absc.sort() ;
  absc = absc[(absc >=0) & (absc <= 1)] ;
  
  if (absc[0] != 0) {
    absc.push_front(0) ;
  }
  if (absc[absc.size() - 1] != 1) {
    absc.push_back(1) ;
  }
  for (int i=0; i < absc.size() ; i++) {
    f.push_back(fdist(absc[i],parms,dist)) ;
  }
  
  fmxind = which_max(f) ;
  
  if (f[fmxind] > fdist(arg1,parms,dist) && f[fmxind] > fdist(arg2,parms,dist)) {
    maxind = 0 ;
  } else if (f[fmxind] < fdist(arg2,parms,dist)) {
    maxind = 1 ;
  } else if (fdist(arg1,parms,dist) > f[fmxind]) {
    maxind = 2 ;
  }
  
  //* calculating intersection x for middle section (containing max); using Rcpp function intsct2
  if  (maxind==2 ) {
    if (fmxind < 3) {
      while (fmxind < 3 ) {
        minp = absc[0];
        pintv = absc[1] - minp ;
        abscinc = pintv/2 ;
        newabsc = minp + abscinc ;
        absc.insert(1, newabsc) ;
        f.insert(1,fdist(newabsc,parms,dist)) ;
        fmxind = which_max(f) ;
      }
    } else if (fmxind>=(absc.size() - 3)) {
      while (fmxind>=(absc.size() - 3)) {
        pintv =  1 - absc[absc.size() - 1];
        abscinc = pintv/2 ;
        newabsc = absc[absc.size() - 2] + abscinc ;
        absc.insert(absc.size() - 1, newabsc) ;
        f.insert(absc.size() - 1, fdist(newabsc,parms,dist)) ;
        fmxind = which_max(f) ;
      }
    } 
  } else {
    
    if (fmxind < 2) {
      while (fmxind < 2) {
        minp = absc[0] ;
        pintv = absc[1] - minp ;
        abscinc = pintv/2 ;
        newabsc = minp + abscinc ;
        absc.insert(1, newabsc) ;
        f.insert(1,fdist(newabsc,parms,dist)) ;
        fmxind = which_max(f) ;
      }
    } else if (fmxind > (absc.size() - 4)) {
      while (fmxind > (absc.size() - 4) ) {
        pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
        abscinc = pintv/2 ;
        newabsc = absc[absc.size() - 2] + abscinc ;
        absc.insert(absc.size() - 1, newabsc) ;
        f.insert(absc.size() - 1,fdist(newabsc,parms,dist)) ;
        fmxind = which_max(f) ;
      }
    }
  }
  
  NumericVector absc2 ;
  for (int i=0 ; i < absc.size() ; i++){
    if (f[i] > -50 || i==0 || i==(absc.size() - 1)) {
      absc2.push_back(absc[i]);
    } 
  }
  NumericVector f2 ;
  for (int i=0 ; i < absc.size() ; i++){
    if (f[i] > -50 || i==0 || i==(absc.size() - 1)) {
      f2.push_back(f[i]);
    } 
  }
  DataFrame df = DataFrame::create( Named("absc")=absc2, Named("f")=f2 ) ;
  return df;
}
