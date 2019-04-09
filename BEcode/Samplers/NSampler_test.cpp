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

double fnbin(double n, NumericVector parms2){
  double x = parms2[0], p = parms2[1] ;
  double f = R::dbinom(x,n,p,false) ;
  return f ;
} 


double fxbin(double x, NumericVector parms2){
  double n = parms2[0], p = parms2[1] ;
  double f = R::dbinom(x,n,p,false) ;
  return f ;
} 

double fxpoi(double x, NumericVector parms2){
  double r = parms2[0] ;
  double f = R::dpois(x,r,false) ;
  return f ;
} 
// parms2 = (lambda,Nhosp,p)
double fxpoinbin(double n, NumericVector parms2){
  double r = parms2[0], p = parms2[2] ;
  int x = parms2[1] ;
  double f = R::dpois(x,r,false) * R::dbinom(x,n,p,false) ;
  return f ;
} 

double fndist(double pr, NumericVector parms2, std::string dist){
  double fv = 0 ;
  if (dist=="binomn") {
    fv = fnbin(pr,parms2) ;
  } else if (dist=="binomx") {
    fv = fxbin(pr,parms2) ;
  } else if (dist=="Poi") {
    fv = fxpoi(pr,parms2) ;
  }else if (dist=="Poibinomn") {
    fv = fxpoinbin(pr,parms2) ;
  }
  return fv ;
} 

DataFrame xbinpdf(NumericVector parms2){
  double p = parms2[1] ;
  int n = parms2[0], x0 = round(n * p) ;
  
  IntegerVector xvec = IntegerVector::create(x0) ;
  NumericVector fvec = NumericVector::create(fxbin(x0,parms2)) ;
  
  DataFrame df ;
  
  int xtest = x0 + 1 ;
  double ftest = fxbin(xtest,parms2) ;
  
  while (ftest > 1e-10) {
    xvec.push_back(xtest) ;
    fvec.push_back(ftest) ;
    xtest++ ;
    ftest = fxbin(xtest,parms2) ;
  }
  
  xtest = x0 - 1 ;
  ftest = fxbin(xtest,parms2) ;
  
  while (ftest > 1e-10 && xtest <= n) {
    xvec.push_front(xtest);
    fvec.push_front(ftest) ;
    xtest-- ;
    ftest = fxbin(xtest,parms2) ;
  }
  
  int xsz = xvec.size();
  double ftotal = 0;
  for(int i = 0 ; i < xsz ; ++i) {
    ftotal += fvec[i];
  }
  
  for(int i = 0 ; i < xsz ; ++i) {
    fvec[i] /= ftotal ;
  }
  
  df = DataFrame::create( Named("x")=xvec, Named("p")=fvec ) ;
  return df ;
} 

DataFrame nbinpdf(NumericVector parms2){
  double p = parms2[1] ;
  int x = parms2[0], n0 = round(x / p) ;
  
  IntegerVector nvec = IntegerVector::create(n0) ;
  NumericVector fvec = NumericVector::create(fnbin(n0,parms2)) ;
  
  DataFrame df ;
  
  int ntest = n0 + 1 ;
  double ftest = fnbin(ntest,parms2) ;
  
  while (ftest > 1e-10) {
    nvec.push_back(ntest) ;
    fvec.push_back(ftest) ;
    ntest++ ;
    ftest = fnbin(ntest,parms2) ;
  }
  
  ntest = n0 - 1 ;
  ftest = fnbin(ntest,parms2) ;
  
  while (ftest > 1e-10 && ntest >= x ) {
    nvec.push_front(ntest);
    fvec.push_front(ftest) ;
    ntest-- ;
    ftest = fnbin(ntest,parms2) ;
  }
  
  int nsz = nvec.size();
  double ftotal = 0;
  for(int i = 0 ; i < nsz ; ++i) {
    ftotal += fvec[i];
  }
  
  for(int i = 0 ; i < nsz ; ++i) {
    fvec[i] /= ftotal ;
  }
  df = DataFrame::create( Named("n")=nvec, Named("p")=fvec ) ;
  return df ;
} 

DataFrame xpoipdf(NumericVector parms2){
  double r = parms2[0] ;
  int x0 = round(r) ;
  
  NumericVector fvec(0) ;
  IntegerVector xvec(0) ;
  
  DataFrame df ;
  
  int xtest = x0 + 1 ;
  double ftest = fxpoi(xtest,parms2) ;
  
  while (ftest > 1e-10) {
    xvec.push_back(xtest) ;
    fvec.push_back(ftest) ;
    xtest++ ;
    ftest = fxpoi(xtest,parms2) ;
  }
  
  if (x0 > 0) {
    xtest = x0 - 1 ;
    ftest = fxpoi(xtest,parms2) ;
    
    while (ftest > 1e-10 && xtest >= 0) {
      xvec.push_front(xtest);
      fvec.push_front(ftest) ;
      xtest-- ;
      ftest = fxpoi(xtest,parms2) ;
    }
  }
  int xsz = xvec.size();
  double ftotal = 0;
  for(int i = 0 ; i < xsz ; ++i) {
    ftotal += fvec[i];
  }
  
  for(int i = 0 ; i < xsz ; ++i) {
    fvec[i] /= ftotal ;
  }
  
  df = DataFrame::create( Named("x")=xvec, Named("p")=fvec ) ;
  return df ;
} 

DataFrame xpoinbinompdf(NumericVector parms2){
  double r = parms2[0] ;
  int x0 = round(r) ;
  
  NumericVector fvec(0) ;
  IntegerVector xvec(0) ;
  
  DataFrame df ;
  
  int xtest = x0 + 1 ;
  double ftest = fxpoi(xtest,parms2) ;
  
  while (ftest > 1e-10) {
    xvec.push_back(xtest) ;
    fvec.push_back(ftest) ;
    xtest++ ;
    ftest = fxpoi(xtest,parms2) ;
  }
  
  if (x0 > 0) {
    xtest = x0 - 1 ;
    ftest = fxpoi(xtest,parms2) ;
    
    while (ftest > 1e-10 && xtest >= 0) {
      xvec.push_front(xtest);
      fvec.push_front(ftest) ;
      xtest-- ;
      ftest = fxpoi(xtest,parms2) ;
    }
  }
  int xsz = xvec.size();
  double ftotal = 0;
  for(int i = 0 ; i < xsz ; ++i) {
    ftotal += fvec[i];
  }
  
  for(int i = 0 ; i < xsz ; ++i) {
    fvec[i] /= ftotal ;
  }
  
  df = DataFrame::create( Named("x")=xvec, Named("p")=fvec ) ;
  return df ;
} 

DataFrame npdf(NumericVector parms2, std::string dist){
  DataFrame df ;
    if (dist=="binomn") {
    df = nbinpdf(parms2) ;
  } else if (dist=="binomx") {
    df = xbinpdf(parms2) ;
  } else if (dist=="Poi") {
    df = xpoipdf(parms2) ;
  }
  return df ;
} 

// [[Rcpp::export]]
IntegerVector nSampler(NumericVector & parms2,std::string dist, int nsims){
  
  DataFrame df = npdf(parms2, dist) ;
  IntegerVector nvec = df[0], outvec ;
  NumericVector fvec = df[1] ;
  double uran, fsum ;
  int outn ;
  
  int selint ;
  
  while (outvec.size() < nsims ) {
    uran =  R::runif(0,1) ;
    
    fsum = 0 ;
    selint = 1; 
    int i = 0 ;
    while (selint == 1) {
      fsum += fvec[i] ;
      if (fsum >=  uran) {
        selint = 0 ;
      }
      i++ ;
    }
    
    outn = nvec[i - 2] ;
    outvec.push_back(outn) ;
  }
  return outvec ;
}
/*** R
parms <- c(100,)
  dist <- "Poi"
dist <- "binom"
absc <- c(.05,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9)
  
  start_time <- Sys.time()
  sls <- samples(absc,parms,dist,10000,.75)
  end_time <- Sys.time()
  end_time - start_time
  hist(sls,100,xlim = c(0,.3),freq = F)
  lines(ls,yls*100,col = 'red')
  */
