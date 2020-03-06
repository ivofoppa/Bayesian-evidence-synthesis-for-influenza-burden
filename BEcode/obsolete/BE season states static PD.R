#########################################################################################
### Estimating hospitalization rates/numbers by age group and state, summing over  ######
### states. Detection probability constant (not estimated). This is to investigate ######
### the need of estimation of state level parameter for national estimates         ######
### Ivo Foppa, 6/2019                                                              ######
#########################################################################################
#########################################################################################
library(R2jags)
#########################################################################################
#########################################################################################
bfolder <- 'C:/Users/VOR1/documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

### Load data 
infname <- 'FluSURV-NET-states.RData'
setwd(paste0(bfolder,'BEdata'))
load(infname)

ptname <- "20180719_Detection_probs_bysite_SAS.csv"
ptdata <- read.csv(ptname, header = T)

ptname2 <- "20190320_Detection_probs_for Ivo_v2.csv"
ptdata2 <- read.csv(ptname2, header = T)

statevec <- unique(as.vector(ptdata$state))
#########################################################################################
###  Analysis with constant detection probs (provided by Carrie/Alyssa) #################
#########################################################################################
codaList <- list()

for (agcat in 1:5) {
  codatotarr <- array(0,dim = c(niter*3/5,1))
  for (state in statevec) {
    selind1 <- which(ptdata$agegrp5==agcat & ptdata$state==state)
    pt <- ptdata[selind1,]$prob_d
    
    selind4 <- which(FSNdata$agecat==agcat & FSNdata$state==state & FSNdata$season==5)
    FSNfluhosp <- sum(FSNdata[selind4,]$freq)
    
    #########################################################################################
    ###  Note: These calcuations for outcome P&I   ##########################################
    #########################################################################################
    
    data <- list('FSNfluhosp'=FSNfluhosp,'pt'=pt)
    
    ## Add one to ensure non-zero denominator
    rfluhospinit <- FSNfluhosp /pt
    
    inits <- function(){
      list(
        rfluhosp = rfluhospinit
      )}
    
    variables <- c('fluhosp')
    
    nadapt <- 10000
    niter <- 10000
    
    model.file <- 'BE state check simple.txt'
    
    setwd(paste0(bfolder,'BEmodels'))
    
    j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
    j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
    codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
    
    codatotarr <- codatotarr + codals
  }
  codaList[[agcat]] <- codatotarr
}

#########################################################################################
#########################################################################################
season <- 5 ## 2014/15
estarr <- array(0,dim = c(0,3))
for (agcat in 1:5) {
  popselind <- which(FSNpopdata$season==season & FSNpopdata$agecat==agcat)
  FSNpop <- FSNpopdata$FSNpop[popselind] 
  USpop <- FSNpopdata$USpop[popselind] 
  fluhospls <- as.vector(codaList[[agcat]])/FSNpop*USpop
  est <- round(as.vector(quantile(fluhospls,probs = c(.5,.025,.975))))
  estarr <- rbind(estarr,est,deparse.level = 0)
}

est <- colSums(estarr)
estarr <- rbind(estarr,est,deparse.level = 0)

setwd(paste0(bfolder,'BEresults'))
write.csv(estarr,file = 'US flu hosp from states, 2014-15.csv',col.names = c('AG','Median','2.5th%tile','97.5th%tile'),
          row.names = c('0-4','5-17','18-49','50-64','65+','All'))
#########################################################################################
###  Analysis with constant detection probs (provided by Carrie/Alyssa) #################
###  The same as before, but aggregated over states (National analysis) #################
#########################################################################################
codaList <- list()

for (agcat in 1:5) {
  selind1 <- which(ptdata2$agegrp5==agcat & ptdata2$season==1415)
  pt <- ptdata2[selind1,]$prob_d
  
  selind4 <- which(FSNdata$agecat==agcat & FSNdata$season==5)
  FSNfluhosp <- sum(FSNdata[selind4,]$freq)
  
  #########################################################################################
  ###  Note: These calcuations for outcome P&I   ##########################################
  #########################################################################################
  
  data <- list('FSNfluhosp'=FSNfluhosp,'pt'=pt)
  
  ## Add one to ensure non-zero denominator
  rfluhospinit <- FSNfluhosp /pt
  
  inits <- function(){
    list(
      rfluhosp = rfluhospinit
    )}
  
  variables <- c('fluhosp')
  
  nadapt <- 10000
  niter <- 10000
  
  model.file <- 'BE state check simple.txt'
  
  setwd(paste0(bfolder,'BEmodels'))
  
  j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  codaList[[agcat]] <- codals
}

#########################################################################################
#########################################################################################
season <- 5 ## 2014/15
estarr <- array(0,dim = c(0,3))
for (agcat in 1:5) {
  popselind <- which(FSNpopdata$season==season & FSNpopdata$agecat==agcat)
  FSNpop <- FSNpopdata$FSNpop[popselind] 
  USpop <- FSNpopdata$USpop[popselind] 
  fluhospls <- as.vector(codaList[[agcat]])/FSNpop*USpop
  est <- round(as.vector(quantile(fluhospls,probs = c(.5,.025,.975))))
  estarr <- rbind(estarr,est,deparse.level = 0)
}

est <- colSums(estarr)
estarr <- rbind(estarr,est,deparse.level = 0)

