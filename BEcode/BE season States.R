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
infname <- 'FluSURV-NET-states.RData'
setwd(paste0(bfolder,'BEdata'))
load(infname)
nseas <- 6
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
agcat <- 5
sensdatasel <- list(sensdata[[1]][[agcat]],sensdata[[2]][[agcat]])
agecatls <- agecatlist[[agcat]]

selind2 <- which(FSNpopdata$agecat==agcat)
agpopsel <- FSNpopdata[selind2,]

#########################################################################################
###  Analyses by state ##################################################################
#########################################################################################
statels <- unique(as.vector(FSNtestdata$state))

st <- "CA"

selind1 <- which(FSNtestdata$agecat==agcat & FSNtestdata$state==st)
FSNtestdatasel <- data.frame(FSNtestdata[selind1,])

selind3 <- which(OSHdata$agecat==agcat & OSHdata$state==st)
oshdatsel <- OSHdata[selind3,]

selind4 <- which(FSNdata$agecat==agcat & FSNdata$state==st)
FSNdatasel <- FSNdata[selind4,]

FSNfluhospls <- FSNdatasel$freq[which(FSNdatasel$died==0 & FSNdatasel$season <= nseas)]
FSNfludeathls <- FSNdatasel$freq[which(FSNdatasel$died==1 & FSNdatasel$season <= nseas)]

dls <- unique(FSNtestdatasel$died)

nttypearr <- array(0,dim = c(2,4,nseas))
ntotarr <- array(0,dim = c(2,nseas))
testposarr <- array(0,dim = c(2,3,nseas))
poshls <- NULL
#########################################################################################
###  Note: These calcuations for outcome P&I   ##########################################
#########################################################################################
for (seas in 1:nseas){
  nttype <- array(0,dim = c(2,4))
  ntot <- c(0,0)
  testpos <- array(0,dim = c(2,3))
  for (d in dls){
    selind2 <- which(FSNtestdatasel$died==d & FSNtestdatasel$season==seas)
    selind2b <- which(FSNtestdatasel$outcome==d & FSNtestdatasel$season==seas)
    ds <- FSNtestdatasel[selind2,]
    dsb <- FSNtestdatasel[selind2b,]
    nttype[d+1,1] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==1)]) # PCR
    nttype[d+1,2] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==2)]) # RIDT
    nttype[d+1,3] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==3)]) # Other/unknown
    nttype[d+1,4] <- sum(ds$freq[which(ds$TestedFlu!=1)]) # not tested
    
    ntot[d+1] <- nttype[d+1,1] + nttype[d+1,2] + nttype[d+1,3] + nttype[d+1,4]
    
    testpos[d+1,1] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==1 & ds$TestResult==1)]) # PCR
    testpos[d+1,2] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==2 & ds$TestResult==1)]) # Other/unknown
    testpos[d+1,3] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==3 & ds$TestResult==1)]) # Other/unknown
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
  selind5 <- which(oshdatsel$season==seas)
  ds3 <- oshdatsel[selind5,]
  
  poshls[seas] <- sum(ds3$pi[which(ds3$osh==1)])/(sum(ds3$pi))
}

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
rfluhosplsinit <- FSNfluhospls /(1 - ptarrinit[1,])
rfludeathlsinit <- FSNfludeathls/(1 - poshls)/(1 - ptarrinit[2,])

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
