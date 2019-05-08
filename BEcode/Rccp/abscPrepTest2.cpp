#include <Rcpp.h>
#include <cmath>
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

double lik(double pr, int parmind, NumericVector & parms, DataFrame & data) {
  NumericVector x = data[0], cases=data[1], controls=data[2], lvec ;
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

double fupperhull(double p, DataFrame & df, int parmind, NumericVector & parms, DataFrame & data, 
                  double zval, int maxind, int zind) {
  int selindp ;
  NumericVector absc = df[0], abscz = absc, l = df[1], lz, absczpart ;
  double l0,l1,l2, uh,p0,p1,p2, a ;
  
  abscz.push_back(zval) ;
  std::sort(abscz.begin(),abscz.end());
  
  for (int i = 0; i < abscz.length() ; i++) {
    lz.push_back(lik(abscz[i], parmind, parms, data)) ;
  }
  
  absczpart = abscz[abscz < p] ;
  selindp = absczpart.length() - 1 ;
  
  if ((abscz[selindp + 1]==zval && zind==0) || (abscz[selindp + 1] < zval )){
    uh = lz[selindp + 1] ;
  } else if (abscz[selindp + 1] == zval && zind==1 && maxind!=0){
    l0 = lz[selindp - 1] ;
    l1 = lz[selindp] ;
    p0 = abscz[selindp - 1] ;
    p1 = abscz[selindp] ;
    a = (l1 - l0)/(p1 - p0) ;
    uh = l1 + a*(p - p1) ;
  } else if (abscz[selindp] == zval && zind==0 && maxind!=0){
    l1 = lz[selindp + 1] ;
    l2 = lz[selindp + 2] ;
    p1 = abscz[selindp + 1] ;
    p2 = abscz[selindp + 2] ;
    a = (l2 - l1)/(p2 - p1) ;
    uh = l1 + a*(p - p1) ;
  } else {
    uh = lz[selindp] ;
  }
  return uh ;
}

int maxindC(DataFrame & df, int parmind, NumericVector & parms, DataFrame & data){
  int lmxind, maxind=3 ;
  double arg1, arg2, lmax, abscmax ;
  NumericVector absc = df[0], l = df[1] ;
  lmax = max(l) ;
  lmxind = which_max(l) ;
  abscmax = absc[lmxind] ;
  
  arg1 = abscmax - 1e-9 ;
  arg2 = abscmax + 1e-9 ;
  
  if (lmax > lik(arg1,parmind, parms, data) && lmax > lik(arg2,parmind, parms, data)) {
    maxind = 0 ;
  } else if (lmax < lik(arg2,parmind, parms, data)) {
    maxind = 1 ;
  } else if (lik(arg1,parmind, parms, data) > lmax) {
    maxind = 2 ;
  }
  
  return maxind ;
} 

double intsct(DataFrame & df, int parmind, NumericVector & parms, DataFrame & data, int maxind){
  
  NumericVector zabsc(4), absc = df[0], f = df[1] ;
  int fmxind = which_max(f) ;
  double p = 1,p1,p2,p3,p4,a = 1,b = 0,f1 = 1,f2 = 1,f3 = 1,f4 = 1 ;
  NumericVector out(zabsc.length() - 3) ;
  
  if (maxind == 2) {
    for (int i=0 ; i < 4 ; i++) {
      zabsc[i] = absc[(fmxind - 2 + i)] ;
    }
  } else {
    for (int i=0 ; i < 4 ; i++) {
      zabsc[i] = absc[(fmxind - 1 + i)] ;
    } 
  }
  
  p1 = zabsc[0] ;
  p2 = zabsc[1] ;
  p3 = zabsc[2] ;
  p4 = zabsc[3] ;
  
  f1 = lik(p1,parmind, parms, data);
  f2 = lik(p2,parmind, parms, data);
  f3 = lik(p3,parmind, parms, data);
  f4 = lik(p4,parmind, parms, data);
  
  a = (f2 - f1)/(p2 - p1) ;
  b = (f4 - f3)/(p4 - p3) ;
  
  p = (f3 - f2 - p3*b + p2*a)/(a - b) ;
  
  return p ;
}

int fzind(DataFrame df, int parmind, NumericVector & parms, DataFrame & data, double zval) {
  double arg ;
  int zind = 0 ;
  NumericVector absc = df[0], f = df[1] ;
  
  arg = zval - 1e-9 ;
  
  if (lik(arg,parmind, parms, data) > lik(zval,parmind, parms, data)) {
    zind = 1 ;
  } else {
    zind = 0 ;
  }
  
  return zind ;
}

List abscPrep(NumericVector & absc0, int parmind, NumericVector & parms, DataFrame & data, double rtio) {
  List outList ;
  DataFrame df ;
  NumericVector l, absc = absc0, abscpart ;
  double abscincr = 0, abscnew = 0, zval, lmax, c1, c2 ;
  int maxind, zind, selind ;
  
  absc = unique(absc) ;
  std::sort(absc.begin(),absc.end()) ;
  
  for (int i = 0; i < absc.length() ; i++) {
    l.push_back(lik(absc[i], parmind, parms, data)) ;
  }
  
  int fmxind = which_max(l) ;
  lmax = max(l) ;
  
  while ((fmxind == 0) || (fmxind == 1) || (exp(l[0] - lmax) > 1e-20)) {
    abscincr = (absc[absc.size() - 1] - absc[0])/10 ;
    abscnew = absc[0] - abscincr ;
    absc.push_front(abscnew) ;
    l.push_front(lik(abscnew, parmind, parms, data)) ;
    fmxind = which_max(l) ;
    lmax = max(l) ;
  }
  
  while ((fmxind == l.size() - 1) || (fmxind == l.size() - 2) || (exp(l[l.size() - 1] - lmax) > 1e-20)) {
    abscincr = (absc[absc.size() - 1] - absc[0])/10 ;
    abscnew = absc[absc.size() - 1] + abscincr ;
    absc.push_back(abscnew) ;
    l.push_back(lik(abscnew, parmind, parms, data)) ;
    fmxind = which_max(l) ;
    lmax = max(l) ;
  }
  
  fmxind = which_max(l) ;

  df = DataFrame::create( Named("absc")=absc, Named("l")=l ) ;
  maxind = maxindC(df, parmind, parms, data) ;
  zval = intsct(df,parmind,parms,data,maxind) ;
  zind = fzind(df,parmind,parms,data,zval) ;
  
  c1 = fupperhull(zval + 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
  c2 = fupperhull(zval - 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
  
  while ((std::max(exp(c1 - lmax),exp(c2 - lmax))/exp(lik(zval,parmind,parms,data) - lmax) > rtio)) {
    abscpart = absc[absc <zval] ;
    selind = abscpart.size() ;
    absc.insert(selind, zval) ;
    l.insert(selind,lik(zval,parmind,parms,data)) ;
    
    lmax = max(l) ;
    
    df = DataFrame::create( Named("absc")=absc, Named("l")=l ) ;
    maxind = maxindC(df, parmind, parms, data) ;
    zval = intsct(df,parmind,parms,data,maxind) ;
    zind = fzind(df,parmind,parms,data,zval) ;
    
    c1 = fupperhull(zval + 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
    c2 = fupperhull(zval - 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
  }

  outList = List::create(Named("df") = df, Named("parmind") = parmind, Named("parms") = parms, Named("zval") = zval, Named("maxind") = maxind, Named("zind") = zind) ;

  return outList ;
}

// [[Rcpp::export]]
List abscGroome(List & outList, DataFrame data, int parmind, int size, double rtio) {
  DataFrame df = as<DataFrame>(outList["df"]) ;
  NumericVector parms = outList["parms"], absc = df[0], l = df[1], zabsc, ldiffvec1(0), ldiffvec2(0), abscadd, abscpart ;
  double lmax = max(l), zval = outList["zval"], abscincr, abscnew, ldiff, c1, c2 ;
  int maxind = outList["maxind"], zind = outList["zind"], fmxind = which_max(l), ldiffind1, ldiffind2, selind ;
  
  fmxind = which_max(l) ;
  lmax = max(l) ;
  
  c1 = fupperhull(zval + 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
  c2 = fupperhull(zval - 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
  
  while ((std::max(exp(c1 - lmax),exp(c2 - lmax))/exp(lik(zval,parmind,parms,data) - lmax) > rtio)) {
    abscpart = absc[absc <zval] ;
    selind = abscpart.size() ;
    absc.insert(selind, zval) ;
    l.insert(selind,lik(zval,parmind,parms,data)) ;
    
    lmax = max(l) ;
    
    df = DataFrame::create( Named("absc")=absc, Named("l")=l ) ;
    maxind = maxindC(df, parmind, parms, data) ;
    zval = intsct(df,parmind,parms,data,maxind) ;
    zind = fzind(df,parmind,parms,data,zval) ;
    
    c1 = fupperhull(zval + 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
    c2 = fupperhull(zval - 1e-6,df, parmind, parms, data, zval, maxind, zind) ;
  }
  
  fmxind = which_max(l) ;
  
  while (fmxind == 0) {
    abscincr = absc[1] - absc[0] ;
    abscnew = absc[0] - abscincr ;
    absc.push_front(abscnew) ;
    l.push_front(lik(abscnew, parmind, parms, data)) ;
    fmxind = which_max(l) ;
  }
  
  while (fmxind == 1) {
    abscincr = (absc[1] - absc[0])/2 ;
    abscnew = absc[0] + abscincr ;
    absc.insert(1, abscnew) ;
    l.insert(1,lik(abscnew,parmind,parms,data)) ;
    fmxind = which_max(l) ;
  }
  
  while (fmxind == (l.size() - 1)) {
    abscincr = absc[absc.size() - 1] - absc[absc.size() - 2] ;
    abscnew = absc[absc.size() - 1] + abscincr ;
    absc.push_back(abscnew) ;
    l.push_back(lik(abscnew, parmind, parms, data)) ;
    fmxind = which_max(l) ;
  }
  
  while (fmxind == (l.size() - 2)) {
    abscincr = (absc[absc.size() - 1] - absc[absc.size() - 2])/2 ;
    abscnew = absc[absc.size() - 2] + abscincr ;
    absc.insert(l.size() - 2, abscnew) ;
    l.insert(l.size() - 2,lik(abscnew,parmind,parms,data)) ;
    fmxind = which_max(l) ;
  }
  
  while(l.size() > size) {
    fmxind = which_max(l) ;
    
    ldiffvec1 = NumericVector::create() ;
    ldiffvec2 = NumericVector::create() ;
    
    for (int k = 1 ; k < fmxind ; k++ ) {
      ldiff = (exp(l[k+1] - lmax) - exp(l[k-1] - lmax))/(absc[k+1] - absc[k-1]) ;
      ldiffvec1.push_back(ldiff) ;
    }
    
    for (int k = fmxind + 1 ; k < (l.size() - 2) ; k++ ) {
      ldiff = abs((exp(l[k+1] - lmax) - exp(l[k-1] - lmax))/(absc[k+1] - absc[k-1])) ;
      ldiffvec2.push_back(ldiff) ;
    }
    
    ldiffind1 = which_min(ldiffvec1) + 1 ;
    ldiffind2 = which_min(ldiffvec2) + fmxind + 1 ;
    
    if (min(ldiffvec1) < min(ldiffvec2)) {
      absc.erase(ldiffind1) ;
      l.erase(ldiffind1) ;
    } else {
      absc.erase(ldiffind2) ;
      l.erase(ldiffind2) ;
    }
  }
  
  df = DataFrame::create( Named("absc")=absc, Named("l")=l ) ;
  maxind = maxindC(df, parmind, parms, data) ;
  zval = intsct(df,parmind,parms,data,maxind) ;
  zind = fzind(df,parmind,parms,data,zval) ;
  List TotList = List::create(Named("df") = df, Named("parmind") = parmind, Named("parms") = parms, Named("zval") = zval, Named("maxind") = maxind, Named("zind") = zind) ;
  
  return TotList ;
}

// [[Rcpp::export]]
List TotPrep(NumericVector & absc, int parmind, NumericVector & parms, DataFrame & data, double rtio) {
  List TotList ;
  
  TotList = abscPrep(absc, parmind, parms, data, rtio) ;
  
  return TotList ;
}

