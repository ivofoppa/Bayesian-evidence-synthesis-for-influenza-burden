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

double lik(double pr,NumericVector & parms,int parmind,DataFrame & data) {
  NumericVector x = data[0],cases=data[1],controls=data[2],lvec ;
  double lv ;
  if (parmind==1) {
    for (int k = 0 ; k < cases.size() ; k++) {
      lv = cases[k]*(pr + parms[1] * x[k] - log(1 + exp(pr + parms[1] * x[k]))) + controls[k]*(- log(1 + exp(pr + parms[1] * x[k])))  ;
      lvec.push_back(lv) ;
    }
  } else if (parmind==2) {
    for (int k = 0 ; k < cases.size() ; k++) {
      lv = cases[k]*(parms[0] + pr * x[k] - log(1 + exp(parms[0] + pr * x[k]))) + controls[k]*(- log(1 + exp(parms[0] + pr * x[k])))  ;
      lvec.push_back(lv) ;
    }
  }
  return sum(lvec) ;
}

// [[Rcpp::export]]
NumericVector fliksum(DataFrame & df,NumericVector & parms,int parmind,DataFrame & data,double zval,int maxind,int zind) {
  NumericVector absc = df[0],l,lnorm,lvec(0),pvec(0) ;
  double l0,l1,p0,p1,a,lmax ;
  
  absc.push_back(zval) ;
  std::sort(absc.begin(),absc.end());
  
  for (int i = 0; i < absc.length() ; i++) {
    l.push_back(lik(absc[i],parms,parmind,data)) ;
  }
  
  lmax = max(l) ;
  
  for (int i = 0 ; i < l.length() ; i++) {
    lnorm.push_back(l[i] - lmax) ;
  }
  
  for (int k = 0 ; k < (absc.size() - 1) ; k++) {
    
    if ((absc[k + 1]==zval && zind==0) || (absc[k + 1] < zval )) {
      lvec.push_back(exp(lnorm[k + 1])*(absc[k + 1] - absc[k])) ;
    } else if (absc[k + 1] == zval && zind==1 && maxind!=0){
      l0 = lnorm[k - 1] ;
      l1 = lnorm[k] ;
      p0 = absc[k - 1] ;
      p1 = absc[k] ;
      a = (l1 - l0)/(p1 - p0) ;
      lvec.push_back(exp(l1)*(exp(a*(zval - p1)) - 1)/a) ;
    } else if (absc[k] == zval && zind==0 && maxind!=0) {
      l0 = lnorm[k + 1] ;
      l1 = lnorm[k + 2] ;
      p0 = absc[k + 1] ;
      p1 = absc[k + 2] ;
      a = (l1 - l0)/(p1 - p0) ;
      lvec.push_back(exp(l0)*(1 - exp(a*(zval - p0)))/a) ;
    } else {
      lvec.push_back(exp(lnorm[k])*(absc[k + 1] - absc[k])) ;
    }
  }

  return lvec ;
}
