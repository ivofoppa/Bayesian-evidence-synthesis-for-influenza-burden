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
double fupperhull(double p,DataFrame & df,NumericVector & parms,int parmind,DataFrame & data,double zval,int maxind,int zind) {
  int selindp ;
  NumericVector absc = df[0],l,abscpart ;
  double l0,l1,l2,uh,p0,p1,p2,a ;
  
  absc.push_back(zval) ;
  std::sort(absc.begin(),absc.end());
  
  for (int i = 0; i < absc.length() ; i++) {
    l.push_back(lik(absc[i],parms,parmind,data)) ;
  }
  
  abscpart = absc[absc < p] ;
  selindp = abscpart.length() - 1 ;
  
  if ((absc[selindp + 1]==zval && zind==0) || (absc[selindp + 1] < zval )){
    uh = l[selindp + 1] ;
  } else if (absc[selindp + 1] == zval && zind==1 && maxind!=0){
    l0 = l[selindp - 1] ;
    l1 = l[selindp] ;
    p0 = absc[selindp - 1] ;
    p1 = absc[selindp] ;
    a = (l1 - l0)/(p1 - p0) ;
    uh = l1 + a*(p - p1) ;
  } else if (absc[selindp] == zval && zind==0 && maxind!=0){
    l1 = l[selindp + 1] ;
    l2 = l[selindp + 2] ;
    p1 = absc[selindp + 1] ;
    p2 = absc[selindp + 2] ;
    a = (l2 - l1)/(p2 - p1) ;
    uh = l1 + a*(p - p1) ;
  } else {
    uh = l[selindp] ;
  }
  return uh ;
}