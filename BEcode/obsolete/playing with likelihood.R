N <- 1000
x <- 479

lik <- function(p){
  choose(N,x)*p^x*(1-p)^(N-x)
}

likls <- sapply(seq(.4,.55,.001), lik)
plot(likls,type = 'l')
plot(log(likls),type = 'l')
