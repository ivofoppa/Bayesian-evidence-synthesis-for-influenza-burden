nsim <- 100000
pls <- NULL


while (length(pls) < nsim) {
  rn1 <- runif(1,.05,.4)
  rn2 <- runif(1)
  if (rn2 < exp(fbin(rn1))/.1) {
    rn1 <- rn1*dbeta(rn1,.5,.5)
    pls <- c(pls,rn1)
  }
}

mean(pls)

nsim2 <- 1e+8
x1 <- rbinom(nsim,n,.2)

p1 <- x1/n

mean(p1)
