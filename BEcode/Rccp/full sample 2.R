library(Rcpp)

setwd('C:/Users/vor1/Dropbox/Misc work/Bayesian debunk/BD code')

sourceCpp('lik.cpp')
sourceCpp('abscPrep.cpp')
sourceCpp('maxindC.cpp')
sourceCpp('intsct.cpp')
sourceCpp('fzind.cpp')

sourceCpp('flowerhull.cpp')
sourceCpp('fupperhull.cpp')

sourceCpp('fliksum.cpp')
sourceCpp('fpsample.cpp')
### This code is for one age group, only for season 2014-15.
### The age group has to be selected; separately for the two data sets!
N1base <- 1
N0base <- 2

Nsize <- 100000

N1 <- round(N1base*Nsize)
N0 <- round(N0base*Nsize)

Nls <- c(N0,N1)

phi <- 2;
ccratio <- 2

mu0 <- 1000
lambda0 <- mu0/N0

nu1 <- mu0*ccratio
gamma <- nu1/N1 ## Sampling ratio

b1 <- log(phi)

for (exp in c(0,1)) {
  eval(parse(text = paste0('N <- N',exp)))
  lambda <- lambda0*exp(b1*exp)
  mu <- lambda*N
  assign(paste0('lambda',exp),lambda)
  assign(paste0('mu',exp),mu)
  assign(paste0('nu',exp),gamma*N)
}
#########################################################################################
### Ensuring the right control-case ratio  
#########################################################################################
nu <- sapply(c(0,1), function(exp) eval(parse(text = paste0('nu',exp))))
mu <- sapply(c(0,1), function(exp) eval(parse(text = paste0('mu',exp))))

rtio <- sum(nu)/sum(mu)
rtio2 <- ccratio/rtio

for (exp in c(0,1)){
  eval(parse(text = paste0('nu',exp,' <- nu[exp + 1]*rtio2')))
}

absc0 <- c(-.4,-.5,-.1,.4,.6)
parms0 <- c(-.7,.5)

#########################################################################################
#########################################################################################
nsims <- 10000
### Simulated case numbers
casearr <- array(0,dim = c(2,nsims))
for (exp in c(0,1)){
  eval(parse(text = paste0('mu <- mu',exp)))
  row <- rpois(nsims,mu)
  casearr[exp + 1,] <- row
}
### Simulated control numbers
controlarr <- array(0,dim = c(2,nsims))
for (exp in c(0,1)){
  eval(parse(text = paste0('nu <- nu',exp)))
  row <- rpois(nsims,nu)
  controlarr[exp + 1,] <- row
}

k <- 1
casels <- casearr[,k]
controlls <- controlarr[,k]

data <- data.frame(exp = c(0,1), cases = casels, controls = controlls)
parmind <- 2

crit <- .5
parms <- parms0
mcmcarr <- analarr <- array(0,dim = c(0,2))

superList <- list()
for (parmind in 1:2) {
  totList <- TotPrep(absc0,parmind,parms,data,2)
  superList[[parmind]] <- totList
}

while (length(mcmcarr[,1]) < nsims) {
  k <- length(mcmcarr[,1]) + 1
  casels <- casearr[,k]
  controlls <- controlarr[,k]
  
  data <- data.frame(exp = c(0,1), cases = casels, controls = controlls)
  # dataast <- data
  # dataast$cases[which(dataast$exp==0)] <- dataast$cases[which(dataast$exp==0)] + 1
  # dataast$controls[which(dataast$exp==1)] <- dataast$controls[which(dataast$exp==1)] + 1
  # 
  # data <- dataast
  # 
  mod1 <- glm(cbind(cases,controls) ~ exp,family = binomial(), data = data)
  summod1 <- summary(mod1)
  ests <- as.vector(summod1$coefficients[,1])
  
  analarr <- rbind(analarr,ests, deparse.level = 0)
  
  
  for (parmind in 1:2) {
    totList <- superList[[parmind]]
    totList <- abscGroome(superlist[[parmind]],data,parmind,15,2)
    superList[[parmind]] <- totList
    
    dfr <- as.data.frame(totList['df'])
    absc <- dfr[,1]
    l <- dfr[,2]
    
    mx <- max(l)
    zval <- unlist(totList['zval'])
    maxind <- unlist(totList['maxind'])
    zind <- unlist(totList['zind'])
    lvec <- fliksum(dfr,parmind,parms,data,zval,maxind,zind)
    
    acc <- 0
    while (acc==0) {
      uran <- runif(1)
      psmp <- fpsample(dfr,parmind,parms,data,zval,maxind,zind,lvec)
      
      lh <- exp(flowerhull(psmp,dfr) - mx)
      uh <- exp(fupperhull(psmp,dfr, parmind, parms, data, zval, maxind, zind) - mx)
      
      if (uran < lh/uh) {
        acc <- 1
      } else {
        lval <- exp(lik(psmp,parmind,parms,data) - mx)
        if (uran < lval/uh) {
          acc <- 1
        } else if (lval/uh < crit) {
          absc1 <- sort(unique(c(psmp,absc)))
          l1 <- sapply(absc1,function(x) )
          
          totList["absc"] <- 
          totList <- abscGroome(totList,data,parmind,15,2)
          
          dfr <- as.data.frame(totList['df'])
          absc <- dfr[,1]
          l <- dfr[,2]
          
          mx <- max(l)
          zval <- unlist(totList['zval'])
          maxind <- unlist(totList['maxind'])
          zind <- unlist(totList['zind'])
          lvec <- fliksum(dfr,parmind,parms,data,zval,maxind,zind)
        } 
      }
    }
    parms[parmind] <- psmp
  }
  mcmcarr <- rbind(mcmcarr,parms, deparse.level = 0)
  if (length(mcmcarr[,1])%%10000==0) {
    cat(length(mcmcarr[,1]),' iterations!\n')
  }
}

b0ls <- mcmcarr[,1]
b1ls <- mcmcarr[,2]
b0glmls <- analarr[,1]
b1glmls <- analarr[,2]

cat("mcmc b0: ", round(mean(b0ls),digits = 4), "  OR: ", round(exp(mean(b1ls)),digits =4))
cat("glm b0: ", round(mean(b0glmls),digits = 4), "  OR: ", round(exp(mean(b1glmls)),digits =4))

mcmcarra2 <- mcmcarr
analarra2 <- analarr


# absc <- absclist[[parmind]]
parmind <- 2
parms <- parms0
totList <- TotPrep(absc0,parmind,parms,data,2)

totList <- abscGroome(totList,data,parmind,15,2)

df <- as.data.frame(totList['df'])
# absc <- unique(df[,1])
# l <- unique(df[,2])
# df <- data.frame(absc,l)

mx <- max(l)
zval <- unlist(totList['zval'])
maxind <- unlist(totList['maxind'])
zind <- unlist(totList['zind'])

xls <- seq(df[2,1],tail(df[,1],2)[1],.001)

lls <- sapply(xls, function(x) lik(x,parmind,parms,data))
lhls <- sapply(xls, function(x) flowerhull(x,df))
uhls <- sapply(xls, function(x) fupperhull(x,df, parmind, parms, data, zval, maxind, zind))

plot(xls,lls,type = 'l',ylim = c(min(lhls),max(uhls)))
lines(xls,lhls,col = 'blue')
lines(xls,uhls,col = 'red')

mx <- max(lls)
likls <- exp(lls - mx)

explhls <- exp(lhls - mx)
expuhls <- exp(uhls - mx)

plot(xls,likls,type = 'l',ylim = c(0,max(expuhls)))
lines(xls,explhls,col = 'blue')
lines(xls,expuhls,col = 'red')


