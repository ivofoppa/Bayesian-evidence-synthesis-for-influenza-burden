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
double fpsample(DataFrame df,NumericVector & parms,int parmind,DataFrame & data,double zval,
                int maxind,int zind,NumericVector lvec) {
    int selind,selindz ;
    
    NumericVector absc = df[0],l = df[1],pveccum,pveccumpart,pvec,abscpart ;
    double l0,l1,p0,p1,a,pran,pbase,ptop,prm,pout,lvecsum,maxl = max(l) ;
    
    lvecsum = sum(lvec) ;
    
    abscpart = absc[absc < zval] ;
    selindz = abscpart.size() ;
    absc.insert(selindz,zval) ;
    l.insert(selindz,lik(zval,parms,parmind,data)) ;

    for (int i = 0 ; i < lvec.size() ; i++) {
      pvec.push_back(lvec[i]/lvecsum) ;
      if (i > 0) {
        pveccum.push_back(pveccum[i - 1] + pvec[i]) ;
      } else {
        pveccum.push_back(pvec[i]) ;
      }
    }
    
    pveccum.push_front(0) ;
    
    pran = R::runif(0,1) ;
    
    pveccumpart = pveccum[pveccum < pran] ;
    selind = pveccumpart.length() - 1 ;
    
    pbase = pveccum[selind] ;
    ptop = pveccum[selind + 1] ;
    prm = pran - pbase ;
    
    if (absc[selind + 1] == zval && zind==1 && maxind > 0) {
      l0 = l[selind - 1] - maxl ;
      l1 = l[selind] - maxl ;
      p0 = absc[selind - 1] ;
      p1 = absc[selind] ;
      a = (l1 - l0)/(p1 - p0) ;
      pout = (log(exp(l1 ) + a*prm*sum(lvec)) + a*p1 - l1)/a ;
    } else if (absc[selind] == zval && zind==0 && maxind > 0) {
      l0 = l[selind + 1] - maxl ;
      l1 = l[selind + 2] - maxl ;
      p0 = absc[selind + 1] ;
      p1 = absc[selind + 2] ;
      a = (l1 - l0)/(p1 - p0) ;
      pout = (log(exp(l0 + a*(zval - p0)) + prm*a*sum(lvec)) + a*p0 - l0)/a ;
    } else {
      pout = absc[selind] + (absc[selind + 1] - absc[selind])*prm/(ptop - pbase) ;
    } 
    return pout ;
  }

