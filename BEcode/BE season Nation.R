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
infname <- 'FluSURV-NET_burden.RData'
setwd(paste0(bfolder,'BEdata'))
load(infname)
nseas <- 5
#########################################################################################
###  The following only needed for US pop data; when doing national analysis  ###########
#########################################################################################
# inputdata <- read.csv('Burden inputs All Ages.csv',header = T)
# agcatls <- unique(agls)
# agls <- as.vector(inputdata$ag)
# agls2 <- sapply(agls, function(ag) which(agcatls==ag))
# agls2 <- as.vector(unlist(agls2))
# inputdata$ag <- agls2
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FluSurvdata,sensdata,agseaspop,oshdat
agcat <- 5

selind1 <- which(FluSurvdata$ag==agcat)
FluSurvdatasel <- FluSurvdata[selind1,]

selind2 <- which(agseaspop$ag==agcat)
agseaspopsel <- agseaspop[selind2,]

selind3 <- which(oshdat$ag==agcat)
oshdatsel <- oshdat[selind3,]

sensdatasel <- list(sensdata[[1]][[agcat]],sensdata[[2]][[agcat]])

agecatls <- agecatlist[[agcat]]
  
nttypearr <- array(0,dim = c(2,4,nseas))
ntotarr <- array(0,dim = c(2,nseas))
testposarr <- array(0,dim = c(2,3,nseas))
poshls <- NULL

FSNfluhosptotls <- NULL
FSNfludeathls <- NULL

for (seas in 1:nseas){
  nttype <- array(0,dim = c(2,4))
  ntot <- c(0,0)
  fluhosp <- 0
  fludeath <- 0
  testpos <- array(0,dim = c(2,3))
  for (d in c(0,1)){
    selind2 <- which(FluSurvdatasel$died==d & FluSurvdatasel$season==seas)
    ds <- FluSurvdatasel[selind2,]
    nttype[d+1,1] <- length(which(ds$TestedFlu==1 & ds$TestType==1)) # PCR
    nttype[d+1,2] <- length(which(ds$TestedFlu==1 & ds$TestType==2)) # RIDT
    nttype[d+1,3] <- length(which(ds$TestedFlu==1 & ds$TestType==3)) # Other/unknown
    nttype[d+1,4] <- length(which(ds$TestedFlu!=1)) # not tested
    
    ntot[d+1] <- nttype[d+1,1] + nttype[d+1,2] + nttype[d+1,3] + nttype[d+1,4]
    fluhosp <- fluhosp + length(ds$died)
    fludeath <- fludeath + d*length(ds$died)
    
    testpos[d+1,1] <- length(which(ds$TestedFlu==1 & ds$TestType==1 & ds$TestResult==1)) # PCR
    testpos[d+1,2] <- length(which(ds$TestedFlu==1 & ds$TestType==2 & ds$TestResult==1)) # Other/unknown
    testpos[d+1,3] <- length(which(ds$TestedFlu==1 & ds$TestType==3 & ds$TestResult==1)) # Other/unknown

  }
  FSNfluhosptotls[seas] <- fluhosp
  FSNfludeathls[seas] <- fludeath
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
  selind4 <- which(oshdatsel$season==seas)
  ds <- oshdatsel[selind4,]
  
  poshls[seas] <- sum(ds$freq[which(ds$osh==1)])/(sum(ds$freq))
  
}

FSNfluhospls <- FSNfluhosptotls - FSNfludeathls ## Non-fatal flu hosp.s
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

FSNpopls <- agseaspopsel[1:5,4]

data <- list('FSNfluhospls'=FSNfluhospls,'FSNpopls'=FSNpopls, 'Npop' = USpop,
               'poshls'=poshls,'FSNfludeathls'=FSNfludeathls,
               'nttypearr'=nttypearr,'testposarr'=testposarr,'pcrsens'=pcrsens,
             'lrapidsens'=lrapidsens,'ntotarr'=ntotarr)

ptesttot <- t(sapply(1:2, function(k) sapply(1:5, function(s) sum(nttypearr[k,,s]))))

ptestarr10init <- nttypearr[,1,]/ptesttot
ptestarr20init <- nttypearr[,2,]/ptesttot
ptestarr30init <- nttypearr[,3,]/ptesttot
ptestarr40init <- nttypearr[,4,]/ptesttot

sens1arrinit <- sapply(1:5, function(s) c(pcrsens[1],pcrsens[1]))
logsens2arrinit <- sapply(1:5, function(s) c(lrapidsens[1],lrapidsens[1]))
sens3arrinit <- sapply(1:5, function(s) c(.4,.4))

sensarrinit <- array(0,dim = c(2,3,5))

for (seas in 1:5) {
  for (k in 1:2){
    sensarrinit[k,1,seas] <- sens1arrinit[k,seas]
    sensarrinit[k,2,seas] <- exp(logsens2arrinit[k,seas])
    sensarrinit[k,3,seas] <- sens3arrinit[k,seas]
  }
}

fluposarrinit <- round(testposarr/sensarrinit)

for (k in 1:2) {
  for (t in 1:3) {
    for (seas in 1:5) {
      if (fluposarrinit[k,t,seas] > nttypearr[k,t,seas]) {
        fluposarrinit[k,t,seas] <- nttypearr[k,t,seas]
      }
    }
  }
}


## Add one to ensure non-zero denominator
pfluarrinit <- sapply(1:nseas, function(s) sapply(c(testposarr[,1,s]/sens1arrinit[,s] + testposarr[,2,s]/exp(logsens2arrinit[,s]) + 
                                                      testposarr[,3,s]/sens3arrinit[,s])/rowSums(nttypearr[,1:3,s]+1),
                   function(x) ifelse(x > .5,min(x,.95),max(x,.1))))

rfluhosplsinit <- FSNfluhospls /FSNpopls*3
rfludeathlsinit <- FSNfludeathls/FSNpopls/(1 - poshls )

inits <- function(){
  list(
    rfluhospls = rfluhosplsinit,
    rfludeathls = rfludeathlsinit,
    ptestarr10 = ptestarr10init,
    ptestarr20 = ptestarr20init,
    ptestarr30 = ptestarr30init,
    ptestarr40 = ptestarr40init,
    pfluarr = pfluarrinit,
    fluposarr = fluposarrinit,
    sens1arr = sens1arrinit,
    logsens2arr = logsens2arrinit,
    sens3arr = sens3arrinit
  )}

variables <- c('rfluhospls')
# variables <- c('pt')
  # variables <- c('pt')
  
nadapt <- 10000
niter <- 10000

model.file <- 'BE season.txt'
  
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
