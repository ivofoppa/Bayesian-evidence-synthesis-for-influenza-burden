library(R2jags)
#########################################################################################
#########################################################################################
rm(list = ls())
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
statevec <- unique(as.vector(FSNtestdata$state))
agcat <- 5

sensdatasel <- list(sensdata[[1]][[agcat]],sensdata[[2]][[agcat]])
agecatvec <- agecatlist[[agcat]]

seasoncodaList <- list()

for (seas in 1:nseas) { 
  
  stateredvec <- as.vector(unlist(sapply(statevec, function(st) if (sum(FSNtestdata$freq[which(FSNtestdata$state==st &  FSNtestdata$TestedFlu==1 &
                                                                                             FSNtestdata$died==1 & FSNtestdata$agecat==agcat &
                                                                                             FSNtestdata$season==seas)]) > 0) st)))
  

  selind2 <- which(FSNpopdata$agecat==agcat)
  #########################################################################################
  ###  Analyses by state ##################################################################
  #########################################################################################
  codatotarr <- array(0,dim = c(6000,12))
  
  for (st in stateredvec) {
    selind1 <- which(FSNtestdata$agecat==agcat & FSNtestdata$state==st & FSNtestdata$season==seas)
    FSNtestdatasel <- data.frame(FSNtestdata[selind1,])
    
    selind3 <- which(OSHdata$agecat==agcat & OSHdata$state==st & OSHdata$season==seas)
    oshdatsel <- OSHdata[selind3,]
    
    selind4 <- which(FSNdata$agecat==agcat & FSNdata$state==st & FSNdata$season==seas)
    FSNdatasel <- FSNdata[selind4,]
    
    FSNfluhosp <- FSNdatasel$freq[which(FSNdatasel$died==0)]
    FSNfludeath <- FSNdatasel$freq[which(FSNdatasel$died==1)]
    
    dvec <- unique(FSNtestdatasel$died)
    
    nttypearr <- array(0,dim = c(2,4))
    ntotvec <- array(0,dim = c(2))
    testposarr <- array(0,dim = c(2,3))
    poshvec <- NULL
    #########################################################################################
    ###  Note: These calcuations for outcome P&I   ##########################################
    #########################################################################################
    nttypearr <- array(0,dim = c(2,4))
    ntot <- c(0,0)
    
    testposarr <- array(0,dim = c(2,3))
    for (d in dvec){
        selind2 <- which(FSNtestdatasel$died==d)
        ds <- FSNtestdatasel[selind2,]
        
        nttypearr[d+1,1] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==1)]) # PCR
        nttypearr[d+1,2] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==2)]) # RIDT
        nttypearr[d+1,3] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==3)]) # Other/unknown
        nttypearr[d+1,4] <- sum(ds$freq[which(ds$TestedFlu!=1)]) # not tested
        
        ntot[d+1] <- nttypearr[d+1,1] + nttypearr[d+1,2] + nttypearr[d+1,3] + nttypearr[d+1,4]
        
        testposarr[d+1,1] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==1 & ds$TestResult==1)]) # PCR
        testposarr[d+1,2] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==2 & ds$TestResult==1)]) # Other/unknown
        testposarr[d+1,3] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==3 & ds$TestResult==1)]) # Other/unknown
      }
      #########################################################################################
      ###  nttype according to all hosp/died system ###########################################
      #########################################################################################

      #########################################################################################
      ###  testpos according to all hosp/died system ##########################################
      #########################################################################################
    ### Continue here: the analysis is by season and state, so all data has to be prepared only for one season
    testpos[1,] <- colSums(testpos)
      testposarr[,,seas] <- testpos
      #########################################################################################
      #########################################################################################
      selind5 <- which(oshdatsel$season==seas)
      ds3 <- oshdatsel[selind5,]
      
      poshvec[seas] <- sum(ds3$pi[which(ds3$osh==1)])/(sum(ds3$pi))
    }
    
    FSNfluhosptot <- FSNfluhosp + FSNfludeath ## Non-fatal flu hosp.s
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
    

    data <- list('FSNfluhosp'=FSNfluhosp,
                 'poshvec'=poshvec,'FSNfludeath'=FSNfludeath,
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
    rfluhospinit <- FSNfluhosp /ptarrinit[1,]
    rfludeathinit <- FSNfludeath/(1 - poshls)/ptarrinit[2,]
    
    pfluarrinit <- array(0,dim=c(2,nseas))
    for (seas in 1:6) {
      pfluarrinit[,seas] <- rowSums(fluposarrinit[,,seas])/rowSums(nttypearr[,,seas])
    }
    
    inits <- function(){
      list(
        rfluhosp = rfluhospinit,
        rfludeath = rfludeathinit,
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
    
    # variables <- c('fludeath')
    variables <- c('fluhosp','fludeath')
    # variables <- c('pt')
    # variables <- c('pt')
    
    nadapt <- 10000
    niter <- 10000
    
    model.file <- 'BE season state.txt'
    
    setwd(paste0(bfolder,'BEmodels'))
    
    j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
    j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
    
    # summary(j.samples)
    
    codaarr <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
    codatotarr <- codatotarr + codals
  }
  seasoncodaList[[seas]] <- list(states=stateredvec,codatotls)
}
#########################################################################################
#########################################################################################
setwd(paste0(bfolder,'BEwriteup'))
fname <- paste0('codalsdiffsimpl',agcat,'.RData')
save(codalsdiffsimpl,file = fname)

#########################################################################################
#########################################################################################
#########################################################################################
plot(codalsdiffsimpl$USfluhosp,type = 'l')
