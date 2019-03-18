nsim <- 100000
pls <- NULL

while (length(pls) < nsim) {
  rn1 <- runif(1,.05,.4)
  rn2 <- runif(1)
  if (rn2 < exp(fbin(rn1))/.12) {
    pls <- c(pls,rn1)
  }
}

mean(pls)

rn1 <- .2
x
