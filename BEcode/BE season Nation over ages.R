rm(list = ls())
library(R2jags)
#########################################################################################
#########################################################################################
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

### Load data 
infname <- 'FSNFullData.csv'
setwd(paste0(bfolder,'BEdata'))
dataset <- read.csv(infname)
#########################################################################################
### Reading-in OSH data      ############################################################
#########################################################################################
OSHfname <- "mort2010_17_season_osh.csv"
OSHdata <- read.csv(OSHfname)
#########################################################################################
nseas <- 7

cipcrlist <- list(c(95.0,82,98.7),c(95.0,82,98.7),c(94.1,81.1,98.7),c(94.1,81.1,98.7),c(86.1,79.6,92.7))
cirapidlist <- list(c(66.7,61.3,71.7),c(66.7,61.3,71.7),c(53.9,47.8,59.8),c(53.9,47.8,59.8),c(20.1,8.8,41.4))
sensdata <- list(cipcrlist,cirapidlist)
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
# agcat <- 5
nadapt <- 10000
niter <- 10000

model.file <- 'BE season ag osh.txt'
for (seas in 1:8) {
  for (ag in 1:5){
    selind <- which(dataset$season==seas & dataset$agcat==ag)
    subdataset <- dataset[selind,]
    
    OSHselind <- which(OSHdata$season==seas & OSHdata$agcat==ag)
    subOSHdata <- OSHdata[OSHselind,]
    
    Npideath <- sum(subOSHdata$pi)
    oshpideath <- subOSHdata$pi[which(subOSHdata$osh==1)]
    
    nttypearr <- array(0,dim = c(2,4))
    ntot <- c(0,0)
    testposarr <- array(0,dim = c(2,3))
    for (d in c(0,1)){
      nttypearr[d+1,1] <- subdataset$ntest1[d+1] # PCR
      nttypearr[d+1,2] <- subdataset$ntest2[d+1] # RIDT
      nttypearr[d+1,3] <- subdataset$ntest3[d+1] # Other/unknown
      nttypearr[d+1,4] <- subdataset$ntest0[d+1] # not tested
      
      ntot[d+1] <- nttypearr[d+1,1] + nttypearr[d+1,2] + nttypearr[d+1,3] + nttypearr[d+1,4]
      
      testposarr[d+1,1] <- subdataset$testres1[d+1] # PCR
      testpos[d+1,2] <- subdataset$testres2[d+1] # Other/unknown
      testposarr[d+1,3] <- subdataset$testres3[d+1] # Other/unknown
    }
    FSNfluhosp <- subdataset$Noutcome[1]
    FSNfludeath <- subdataset$Noutcome[2]
    FSNfluhosptot <- FSNfluhosp + FSNfludeath ## Non-fatal flu hosp.s
    #########################################################################################
    #########################################################################################
    #########################################################################################
    sensdatasel <- list(sensdata[[1]][[ag]],sensdata[[2]][[ag]])
    cipcr <- sensdatasel[[1]]/100
    seest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
    pcrsens <- c(cipcr[1],seest)
    #########################################################################################
    #########################################################################################
    #########################################################################################
    cirapid <- sensdatasel[[2]]/100
    
    selograpidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
    lrapidsens <- c(log(cirapid[1]),selograpidest)
    
    FSNpop <- subdataset$FSNpop[1]
    USpop  <- subdataset$USpop[1]
    
    oshpideath <- subOSHdata$pi[which(subOSHdata$osh==1)]
    Npideath <- sum(subOSHdata$pi)
    
    data <- list('FSNfluhosp'=FSNfluhosp,'FSNpop'=FSNpop, 'USpop' = USpop,
                 'oshpideath'=oshpideath,'Npideath'=Npideath,'FSNfludeath'=FSNfludeath,
                 'nttypearr'=nttypearr,'testposarr'=testposarr,'pcrsens'=pcrsens,
                 'lrapidsens'=lrapidsens,'ntot'=ntot)
    
    
  }
}
  #########################################################################################
  ###  Note: These calcuations for outcome P&I   ##########################################
  #########################################################################################
  for (seas in 1:nseas){
    #########################################################################################
    ###  testpos according to all hosp/died system ##########################################
    #########################################################################################

    
  ptesttot <- t(sapply(1:2, function(k) sum(nttypearr[k,])))
  
  ptest10init <- nttypearr[,1]/ptesttot
  ptest20init <- nttypearr[,2]/ptesttot
  ptest30init <- nttypearr[,3]/ptesttot
  ptest40init <- nttypearr[,4]/ptesttot
  
  sens1init <- c(pcrsens[1],pcrsens[1])
  logsens2init <- c(lrapidsens[1],lrapidsens[1])
  sens3init <- c(.4,.4)
  
  sensarrinit <- array(0,dim = c(2,3))
  
    for (k in 1:2){
      sensarrinit[k,1] <- sens1init[k]
      sensarrinit[k,2] <- exp(logsens2init[k])
      sensarrinit[k,3] <- sens3init[k]
    }
  }
  
  fluposarrinit <- round(testposarr/sensarrinit)
  
  for (k in 1:2) {
    for (t in 1:3) {
        if (fluposarrinit[k,t] > nttypearr[k,t]) {
          fluposarrinit[k,t] <- nttypearr[k,t]
        }
      }
  }
  
  
  ## Add one to ensure non-zero denominator
  poshinit <- oshpideath/Npideath
  rfluhospinit <- FSNfluhosp /FSNpop*3
  rfludeathinit <- FSNfludeath/FSNpop/(1 - poshinit)
  
  USfludeathinit <- round(rfludeathinit*USpop)
  USfluhospinit <- round(rfluhospinit*USpop)
  
  ptot <- NULL
  
  ptest10 <- ptest10init
  ptest20 <- ptest20init
  ptest30 <- ptest30init
  ptest40 <- ptest40init
  
  ptestarr <- array(0,dim = c(2,4))
  sensarr <- array(0,dim = c(2,3))
  
  ptinit <- NULL
  for (k in 1:2) {
    ptot[k] <- ptest10[k] + ptest20[k] + ptest30[k] + ptest40[k]
    
    ptestarr[k,1] <- ptest10[k]/ptot[k]
    ptestarr[k,2] <- ptest20[k]/ptot[k]
    ptestarr[k,3] <- ptest30[k]/ptot[k]
    ptestarr[k,4] <- ptest40[k]/ptot[k]
    
    sensarr[k,1] <- sens1init[k]
    sensarr[k,2] <- exp(logsens2init[k])
    sensarr[k,3] <- sens3init[k]
    
    ptinit[k] <- ptestarr[k,1]*sensarr[k,1] + ptestarr[k,2]*sensarr[k,2] + 
      ptestarr[k,3]*sensarr[k,3]
    
  }
  
  fludeathinit <- round(FSNfludeath/ptinit[2])
  fluhospinit <- round(FSNfluhosp/ptinit[1])
  
  inits <- function(){
    list(
      fludeath = fludeathinit,
      fluhosp = fluhospinit,
      rfluhosp = rfluhospinit,
      rfludeath = rfludeathinit,
      ptest10 = ptest10init,
      ptest20 = ptest20init,
      ptest30 = ptest30init,
      ptest40 = ptest40init,
      fluposarr = fluposarrinit,
      sens1 = sens1init,
      logsens2 = logsens2init,
      sens3 = sens3init,
      posh = poshinit
    )}
  
  variables <- c('USfludeath','USfluhosp')
  # variables <- c('USfluhospls')
  # variables <- c('pt')
  # variables <- c('pt')
  
  setwd(paste0(bfolder,'BEmodels'))
  
  j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
  
  # summary(j.samples)
  
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  codalsdiffsimpl <- data.frame(codals)
  #########################################################################################
  #########################################################################################
  setwd(paste0(bfolder,'BEwriteup'))
  fname <- paste0('codalsdiffsimpl',agcat,'.RData')
  save(codalsdiffsimpl,file = fname)
#########################################################################################
#########################################################################################
#########################################################################################
plot(codalsdiffsimpl$USfluhosp,type = 'l')
