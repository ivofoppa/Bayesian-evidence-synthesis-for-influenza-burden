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

double lik(double pr, NumericVector & parms,int parmind,DataFrame & data) {
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

double flowerhull(double p,DataFrame & df) {
  int selind = 0 ;
  NumericVector absc = df[0],l = df[1],abscpart ;
  double l0,l1,p0,p1,lh,xcap = 0,pmax ;
  
  abscpart = absc[absc < p] ;
  pmax = max(abscpart) ;
  
  selind = abscpart.length() - 1 ;
  
  xcap = p - pmax ;
  l0 = l[selind] ;
  l1 = l[selind + 1] ;
  
  p0 = absc[selind] ;
  p1 = absc[selind + 1] ;
  
  lh = l0 + (l1 - l0) * xcap/(p1 - p0) ;
  return lh ;
}

double fupperhull(double p,DataFrame & df, NumericVector & parms,int parmind,DataFrame & data,
                  double zval,int maxind,int zind) {
  int selindp ;
  NumericVector absc = df[0],abscz = absc,l = df[1],lz,absczpart ;
  double l0,l1,l2,uh,p0,p1,p2,a ;
  
  abscz.push_back(zval) ;
  std::sort(abscz.begin(),abscz.end());
  
  for (int i = 0; i < abscz.length() ; i++) {
    lz.push_back(lik(abscz[i],parms,parmind,data)) ;
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

int maxindC(DataFrame & df, NumericVector & parms,int parmind,DataFrame & data){
  int lmxind,maxind=3 ;
  double arg1,arg2,lmax,abscmax ;
  NumericVector absc = df[0],l = df[1] ;
  lmax = max(l) ;
  lmxind = which_max(l) ;
  abscmax = absc[lmxind] ;
  
  arg1 = abscmax - 1e-20 ;
  arg2 = abscmax + 1e-20 ;
  
  if (lmax > lik(arg1,parms,parmind,data) && lmax > lik(arg2,parms,parmind,data)) {
    maxind = 0 ;
  } else if (lmax < lik(arg2,parms,parmind,data)) {
    maxind = 1 ;
  } else if (lik(arg1,parms,parmind,data) > lmax) {
    maxind = 2 ;
  }
  
  return maxind ;
} 

double intsct(DataFrame & df, NumericVector & parms,int parmind,DataFrame & data,int maxind){
  
  NumericVector zabsc(4),absc = df[0],f = df[1] ;
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
  
  f1 = lik(p1,parms,parmind,data);
  f2 = lik(p2,parms,parmind,data);
  f3 = lik(p3,parms,parmind,data);
  f4 = lik(p4,parms,parmind,data);
  
  a = (f2 - f1)/(p2 - p1) ;
  b = (f4 - f3)/(p4 - p3) ;
  
  p = (f3 - f2 - p3*b + p2*a)/(a - b) ;
  
  return p ;
}

int fzind(DataFrame df, NumericVector & parms,int parmind,DataFrame & data,double zval) {
  double arg ;
  int zind = 0 ;
  NumericVector absc = df[0],f = df[1] ;
  
  arg = zval - 1e-20 ;
  
  if (lik(arg,parms,parmind,data) > lik(zval,parms,parmind,data)) {
    zind = 1 ;
  } else {
    zind = 0 ;
  }
  
  return zind ;
}

List abscGroome(List & outList, NumericVector & parms, int parmind, DataFrame data, int size, double rtio) {
  DataFrame df = as<DataFrame>(outList["df"]) ;
  NumericVector absc = df[0], l, zabsc, ldiffvec1(0), ldiffvec2(0), abscadd, abscpart ;
  double lmax = max(l), zval = outList["zval"], abscincr, abscnew, ldiff, c1, c2 ;
  int maxind = outList["maxind"], zind = outList["zind"], fmxind = which_max(l), ldiffind1, ldiffind2, selind ;
  
  for (int k = 0; k < absc.size() ; k++) {
    l.push_back(lik(absc[k],parms,parmind,data)) ;
  }
  
  fmxind = which_max(l) ;
  lmax = max(l) ;
  
  while ((fmxind == 0) || (fmxind == 1) || (exp(l[0] - lmax) > 1e-50)) {
    abscincr = (absc[absc.size() - 1] - absc[0])/absc.size() ;
    abscnew = absc[0] - abscincr ;
    absc.push_front(abscnew) ;
    l.push_front(lik(abscnew, parms, parmind, data)) ;
    fmxind = which_max(l) ;
    lmax = max(l) ;
  }
  
  while ((fmxind == l.size() - 1) || (fmxind == l.size() - 2) || (exp(l[l.size() - 1] - lmax) > 1e-50)) {
    abscincr = (absc[absc.size() - 1] - absc[0])/absc.size() ;
    abscnew = absc[absc.size() - 1] + abscincr ;
    absc.push_back(abscnew) ;
    l.push_back(lik(abscnew, parms, parmind, data)) ;
    fmxind = which_max(l) ;
    lmax = max(l) ;
  }
  
  fmxind = which_max(l) ;
  
  c1 = fupperhull(zval + 1e-20,df, parms, parmind, data, zval, maxind, zind) ;
  c2 = fupperhull(zval - 1e-20,df, parms, parmind, data, zval, maxind, zind) ;
  
  while ((std::max(exp(c1 - lmax),exp(c2 - lmax))/exp(lik(zval,parms,parmind,data) - lmax) > rtio)) {
    abscpart = absc[absc <zval] ;
    selind = abscpart.size() ;
    absc.insert(selind, zval) ;
    l.insert(selind,lik(zval,parms,parmind,data)) ;
    
    lmax = max(l) ;
    
    df = DataFrame::create( Named("absc")=absc, Named("l")=l ) ;
    maxind = maxindC(df, parms, parmind, data) ;
    zval = intsct(df,parms,parmind,data,maxind) ;
    zind = fzind(df,parms,parmind,data,zval) ;
    
    c1 = fupperhull(zval + 1e-20,df, parms, parmind, data, zval, maxind, zind) ;
    c2 = fupperhull(zval - 1e-20,df, parms, parmind, data, zval, maxind, zind) ;
  }
  
  fmxind = which_max(l) ;
  
  while (fmxind == 0) {
    abscincr = absc[1] - absc[0] ;
    abscnew = absc[0] - abscincr ;
    absc.push_front(abscnew) ;
    l.push_front(lik(abscnew, parms, parmind, data)) ;
    fmxind = which_max(l) ;
  }
  
  while (fmxind == 1) {
    abscincr = (absc[1] - absc[0])/2 ;
    abscnew = absc[0] + abscincr ;
    absc.insert(1, abscnew) ;
    l.insert(1,lik(abscnew,parms,parmind,data)) ;
    fmxind = which_max(l) ;
  }
  
  while (fmxind == (l.size() - 1)) {
    abscincr = absc[absc.size() - 1] - absc[absc.size() - 2] ;
    abscnew = absc[absc.size() - 1] + abscincr ;
    absc.push_back(abscnew) ;
    l.push_back(lik(abscnew, parms, parmind, data)) ;
    fmxind = which_max(l) ;
  }
  
  while (fmxind == (l.size() - 2)) {
    abscincr = (absc[absc.size() - 1] - absc[absc.size() - 2])/2 ;
    abscnew = absc[absc.size() - 2] + abscincr ;
    absc.insert(l.size() - 2, abscnew) ;
    l.insert(l.size() - 2,lik(abscnew,parms,parmind,data)) ;
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
  maxind = maxindC(df, parms, parmind, data) ;
  zval = intsct(df,parms,parmind,data,maxind) ;
  zind = fzind(df,parms,parmind,data,zval) ;
  List TotList = List::create(Named("df") = df, Named("parmind") = parmind, Named("zval") = zval, Named("maxind") = maxind, Named("zind") = zind) ;
  
  return TotList ;
}

NumericVector fliksum(DataFrame & df, NumericVector & parms,int parmind,DataFrame & data,double zval,int maxind,int zind) {
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

// [[Rcpp::export]]
DataFrame abscSampler(List & superList,NumericVector & parms,DataFrame & data,int size,double rtio,double crit,int nsims) {
  List inList;
  DataFrame parmSampleArr,df ;
  NumericVector absc,l,lvec,abscpart,parmSample,parmS1,parmS2 ;
  int parmind,maxind,zind,acc,selind ;
  double mx,zval,uran,psmp,lh,uh,lval ; 
  
  for (int sim = 0 ; sim < nsims ; sim++) {
    parmSample = NumericVector::create() ;
    for (parmind = 1 ; parmind <3 ; parmind++) {
      inList = as<List>(superList[parmind - 1]) ;
      inList = abscGroome(inList,parms,parmind,data,size,rtio) ;
      
      df = as<DataFrame>(inList["df"]) ;
      absc = df[0] ;
      l = df[1] ;
      
      mx = max(l) ;
      zval = inList["zval"] ;
      maxind = inList["maxind"] ;
      zind = inList["zind"] ;
      
      acc = 0 ;
      
      while (acc==0) {
       uran = R::runif(0,1) ;
        lvec = fliksum(df,parms,parmind,data,zval,maxind,zind) ;
        psmp = fpsample(df,parms,parmind,data,zval,maxind,zind,lvec) ;
        
        lh = exp(flowerhull(psmp,df) - mx) ;
        uh = exp(fupperhull(psmp,df,parms,parmind,data,zval,maxind,zind) - mx) ;
        
        if (uran < lh/uh) {
          acc = 1 ;
        } else {
          lval = exp(lik(psmp,parms,parmind,data) - mx) ;
          if (uran < lval/uh) {
            acc = 1 ;
          } else if (lval/uh < crit) {
            if (psmp < absc[0]) {
              stop("problemo!") ;
            }
            abscpart = absc[absc < psmp] ;
            selind = abscpart.size() ;
            absc.insert(selind,psmp) ;
            l.insert(selind,lik(psmp,parms,parmind,data)) ;
            mx = max(l) ;
            df = DataFrame::create(Named("absc") = absc,Named("l") = l) ;
            
            maxind = maxindC(df,parms,parmind,data) ;
            zval = intsct(df,parms,parmind,data,maxind) ;
            zind = fzind(df,parms,parmind,data,zval) ;
            
            inList = List::create(Named("df") = df,Named("parmind") = parmind,Named("zval") = zval,Named("maxind") = maxind,Named("zind") = zind) ;
            superList[parmind - 1] = inList ;
          } 
        }
      }
      parms[parmind - 1] = psmp ;
      superList[parmind - 1] = inList ;
      parmSample.push_back(psmp) ;
      if (parmind==1) {
        parmS1.push_back(psmp) ;
      } else {
        parmS2.push_back(psmp) ;
      }
    }
  }
  parmSampleArr = DataFrame::create(Named("b0") = parmS1,Named("b1") = parmS2) ;
  
  return parmSampleArr ;
}


