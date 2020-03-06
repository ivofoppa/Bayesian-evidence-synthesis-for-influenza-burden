alphals <- seq(1,100,.1)
cipcr <- cipcrlist[[5]]/100
mu <- cipcr[1]

betals <- sapply(alphals, function(a) a*(1-mu)/mu)

ciarr <- array(0,dim=c(length(alphals),2))

for (k in seq_along(alphals)) {
  alpha <- alphals[k]
  beta <- betals[k]
  rls <- rbeta(10000,alpha,beta)
  ci <- as.vector(quantile(rls,probs = c(.5,.025,.975)))
  ciarr[k,] <- ci[2:3]
}

sqrerr <- (ciarr[,1]-cipcr[2])^2 + (ciarr[,2]-cipcr[3])^2

selind <- which(sqrerr==min(sqrerr))

alpha <- alphals[selind]
beta <- betals[selind]

rls <- rbeta(100000,alpha,beta)

quantile(rls,probs = c(.5,.025,.975))

cipcr
###################################################################################################
###  "Fitting" a log-Normal instead   #############################################################
###################################################################################################
logseest <- (log(cipcr[3]) - log(cipcr[2]))/2/1.96
logmu <- log(cipcr[1])

logrls <- rnorm(100000,logmu,sd=logseest)
exp(as.vector(quantile(logrls,probs = c(.5,.025,.975))))

cipcr
################## rapid
cirapid <- sensdatasel[[2]]/100

logseest <- (log(cirapid[3]) - log(cirapid[2]))/2/1.96
logmu <- log(cirapid[1])

logrls <- rnorm(100000,logmu,sd=logseest)

exp(as.vector(quantile(logrls,probs = c(.5,.025,.975))))

cirapid
################# beta fit with rapid
alphals <- seq(1,100,.1)

mu <- cirapid[1]

betals <- sapply(alphals, function(a) a*(1-mu)/mu)

ciarr <- array(0,dim=c(length(alphals),2))

for (k in seq_along(alphals)) {
  alpha <- alphals[k]
  beta <- betals[k]
  rls <- rbeta(10000,alpha,beta)
  ci <- as.vector(quantile(rls,probs = c(.5,.025,.975)))
  ciarr[k,] <- ci[2:3]
}

sqrerr <- (ciarr[,1]-cirapid[2])^2 + (ciarr[,2]-cirapid[3])^2

selind <- which(sqrerr==min(sqrerr))

alpha <- alphals[selind]
beta <- betals[selind]

rls <- rbeta(100000,alpha,beta)

quantile(rls,probs = c(.5,.025,.975))

cirapid


