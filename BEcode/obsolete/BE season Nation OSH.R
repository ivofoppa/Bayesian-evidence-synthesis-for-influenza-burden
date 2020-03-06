rm(list = ls())

library(R2jags)
#########################################################################################
#########################################################################################
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

### Load data 
infname <- 'FluSURV-NET-Nation-8seas.RData'
setwd(paste0(bfolder,'BEdata'))
load(infname)
nseas <- 7
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
agcat <- 5

selind1 <- which(Fulldata$agcat==agcat)
Fulldatasel <- Fulldata[selind1,]

selind2 <- which(OSHdata$agcat==agcat)
OSHdatasel <- OSHdata[selind2,]

sensdatasel <- list(sensdata[[1]][[agcat]],sensdata[[2]][[agcat]])

agecatls <- agecatlist[[agcat]]
  
nttypearr <- array(0,dim = c(2,4,nseas))
ntotarr <- array(0,dim = c(2,nseas))
testposarr <- array(0,dim = c(2,3,nseas))
poshls <- NULL

FSNfluhospls <- Fulldatasel$Noutcome[which(Fulldatasel$died==0 & Fulldatasel$season <= nseas)]
FSNfludeathls <- Fulldatasel$Noutcome[which(Fulldatasel$died==1 & Fulldatasel$season <= nseas)]

Npideathls <- oshpideathls <- NULL

#########################################################################################
###  Note: These calcuations for outcome P&I   ##########################################
#########################################################################################
for (seas in 1:nseas){
  nttype <- array(0,dim = c(2,4))
  ntot <- c(0,0)
  testpos <- array(0,dim = c(2,3))
  for (d in c(0,1)){
    selind2 <- which(Fulldatasel$died==d & Fulldatasel$season==seas)
    ds <- Fulldatasel[selind2,]

    nttype[d+1,1] <- ds$ntest1 # PCR
    nttype[d+1,2] <- ds$ntest2 # RIDT
    nttype[d+1,3] <- ds$ntest3 # Other/unknown
    nttype[d+1,4] <- ds$ntest0 # not tested
    
    ntot[d+1] <- nttype[d+1,1] + nttype[d+1,2] + nttype[d+1,3] + nttype[d+1,4]
    
    testpos[d+1,1] <- ds$testres1 # PCR
    testpos[d+1,2] <- ds$testres2 # Other/unknown
    testpos[d+1,3] <- ds$testres3 # Other/unknown
  }
  
  #########################################################################################
  ###  nttype according to all hosp/died system ###########################################
  #########################################################################################
  nttypearr[,,seas] <- nttype
  ntotarr[,seas] <- ntot
  
  #########################################################################################
  ###  testpos according to all hosp/died system ##########################################
  #########################################################################################
  testpos[1,] <- colSums(testpos)
  testposarr[,,seas] <- testpos
  #########################################################################################
  #########################################################################################
  Npideathls[seas] <- sum(OSHdatasel$pi[which(OSHdatasel$season==seas )])
  oshpideathls[seas] <- OSHdatasel$pi[which(OSHdatasel$season==seas & OSHdatasel$osh==1)]
  
}
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

FSNpopls <- sapply(1:nseas, function(s) Fulldatasel$FSNpop[which(Fulldatasel$season==s)][1])
USpopls  <- sapply(1:nseas, function(s) Fulldatasel$USpop[which(Fulldatasel$season==s)][1])

data <- list('FSNfluhospls'=FSNfluhospls,'FSNpopls'=FSNpopls, 'USpopls' = USpopls,
             'Npideathls'=Npideathls,'oshpideathls'=oshpideathls,'FSNfludeathls'=FSNfludeathls,
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
poshlsinit <- oshpideathls/Npideathls

rfluhosplsinit <- FSNfluhospls /FSNpopls*3
rfludeathlsinit <- FSNfludeathls/FSNpopls/(1 - poshlsinit)

USfludeathlsinit <- round(rfludeathlsinit*USpopls)
USfluhosplsinit <- round(rfluhosplsinit*USpopls)


inits <- function(){
  list(
    rfluhospls = rfluhosplsinit,
    rfludeathls = rfludeathlsinit,
    ptestarr10 = ptestarr10init,
    ptestarr20 = ptestarr20init,
    ptestarr30 = ptestarr30init,
    ptestarr40 = ptestarr40init,
    fluposarr = fluposarrinit,
    sens1arr = sens1arrinit,
    logsens2arr = logsens2arrinit,
    sens3arr = sens3arrinit,
    poshls = poshlsinit,
    USfludeathls = USfludeathlsinit,
    USfluhospls = USfluhosplsinit
  )}

variables <- c('USfludeathls','USfluhospls')
# variables <- c('USfluhospls')
# variables <- c('pt')
  # variables <- c('pt')
  
nadapt <- 10000
niter <- 10000

model.file <- 'BE season OSH.txt'
  
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
}
#########################################################################################
#########################################################################################
#########################################################################################
plot(codalsdiffsimpl$USfluhosp,type = 'l')
