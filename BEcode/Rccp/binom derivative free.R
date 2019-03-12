absc <- c(.1,.15,.2,.25,.4,.5)

n <- 100; x <- 20;

f <- function(p){
  dbinom(x,n,p,log = T)
}

# intsct <- function(p1,p2,p3,p4){
#   a <- (f(p2) - f(p1))/(p2 - p1)
#   b <- (f(p4) - f(p3))/(p4 - p3)
#   y <- (f(p3) - f(p1) + p1*a - p3*b)/(a - b)
#   (y - f(p1) + p1*a)/a
# }
# 
# 
intsct <- function(absc0){
  absc <- sort(absc0)
  k <- 1
  outls <- NULL
  while ((k + 3) <= length(absc)){
    
    p1 <- absc[k];p2 <- absc[k + 1];p3 <- absc[k + 2];p4 <- absc[k + 3]
    
    a <- (f(p2) - f(p1))/(p2 - p1)
    b <- (f(p4) - f(p3))/(p4 - p3)
    p <- (f(p3) - f(p1) + p1*a - p3*b)/(a - b)
    outls <- c(outls,p)
    k <- k + 1
  }
  outls
}

intsct(absc)


