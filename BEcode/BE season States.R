FSNfluhosptotls <- FSNfluhospls + FSNfludeathls ## Non-fatal flu hosp.s
#########################################################################################
#########################################################################################
#########################################################################################
cipcr <- sensdatasel[[1]]/100
seest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
pcrsens <- c(cipcr[1],seest)
#########################################################################################
#########################################################################################
#########################################################################################
cirapid <- sensdatasel[[2]]/100

selograpidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
lrapidsens <- c(log(cirapid[1]),selograpidest)

FSNpopls <- agpopsel[1:nseas,4]
USpopls  <- agpopsel[1:nseas,3]

data <- list('FSNfluhospls'=FSNfluhospls,
               'poshls'=poshls,'FSNfludeathls'=FSNfludeathls,
               'nttypearr'=nttypearr,'testposarr'=testposarr,'pcrsens'=pcrsens,
             'lrapidsens'=lrapidsens,'ntotarr'=ntotarr,'nseas'=nseas)

ptesttot <- t(sapply(1:2, function(k) sapply(1:nseas, function(s) sum(nttypearr[k,,s]))))

ptestarr10init <- nttypearr[,1,1:nseas]/ptesttot
ptestarr20init <- nttypearr[,2,1:nseas]/ptesttot
ptestarr30init <- nttypearr[,3,1:nseas]/ptesttot
ptestarr40init <- nttypearr[,4,1:nseas]/ptesttot

sens1arrinit <- sapply(1:nseas, function(s) c(pcrsens[1],pcrsens[1]))
logsens2arrinit <- sapply(1:nseas, function(s) c(lrapidsens[1],lrapidsens[1]))
sens3arrinit <- sapply(1:nseas, function(s) c(.4,.4))

sensarrinit <- array(0,dim = c(2,3,nseas))

for (seas in 1:nseas) {
  for (k in 1:2){
    sensarrinit[k,1,seas] <- sens1arrinit[k,seas]
    sensarrinit[k,2,seas] <- exp(logsens2arrinit[k,seas])
    sensarrinit[k,3,seas] <- sens3arrinit[k,seas]
  }
}

ptestarrinit <- array(0,dim = c(2,4,nseas))

for (seas in 1:nseas) {
  for (k in 1:2){
    ptestarrinit[k,1,seas] <- ptestarr10init[k,seas]
    ptestarrinit[k,2,seas] <- ptestarr20init[k,seas]
    ptestarrinit[k,3,seas] <- ptestarr30init[k,seas]
    ptestarrinit[k,4,seas] <- ptestarr40init[k,seas]
  }
}

ptarrinit <- array(0,dim = c(2,nseas))

for (seas in 1:nseas) {
  ptarrinit[1,seas] <- sum(sensarrinit[1,,seas]*ptestarrinit[1,1:3,seas])
  ptarrinit[2,seas] <- sum(sensarrinit[2,,seas]*ptestarrinit[2,1:3,seas])
}


fluposarrinit <- round(testposarr/sensarrinit)

for (k in 1:2) {
  for (t in 1:3) {
    for (seas in 1:nseas) {
      if (fluposarrinit[k,t,seas] > nttypearr[k,t,seas]) {
        fluposarrinit[k,t,seas] <- nttypearr[k,t,seas]
      }
    }
  }
}


## Add one to ensure non-zero denominator
rfluhosplsinit <- FSNfluhospls /ptarrinit[1,]
rfludeathlsinit <- FSNfludeathls/(1 - poshls)/ptarrinit[2,]

pfluarrinit <- array(0,dim=c(2,nseas))
for (seas in 1:6) {
  pfluarrinit[,seas] <- rowSums(fluposarrinit[,,seas])/rowSums(nttypearr[,,seas])
}

inits <- function(){
  list(
    rfluhospls = rfluhosplsinit,
    rfludeathls = rfludeathlsinit,
    ptestarr10 = ptestarr10init,
    ptestarr20 = ptestarr20init,
    ptestarr30 = ptestarr30init,
    ptestarr40 = ptestarr40init,
    fluposarr = fluposarrinit,
    pfluarr = pfluarrinit,
    sens1arr = sens1arrinit,
    logsens2arr = logsens2arrinit,
    sens3arr = sens3arrinit
  )}

# variables <- c('fludeathls')
variables <- c('fluhospls','fludeathls')
# variables <- c('pt')
  # variables <- c('pt')
  
nadapt <- 10000
niter <- 10000

model.file <- 'BE season state.txt'
  
setwd(paste0(bfolder,'BEmodels'))
  
j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 

summary(j.samples)
  
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  codalsdiffsimpl <- data.frame(codals)
  #########################################################################################
  #########################################################################################
  setwd(paste0(bfolder,'BEwriteup'))
  fname <- paste0('codalsdiffsimpl',agcat,'.RData')
  save(codalsdiffsimpl,file = fname)
}
#########################################################################################
#########################################################################################
#########################################################################################
plot(codalsdiffsimpl$USfluhosp,type = 'l')
