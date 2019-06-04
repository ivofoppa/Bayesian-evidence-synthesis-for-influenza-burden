library(readxl)
#########################################################################################
#########################################################################################
rm(list = ls())
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

seasvec <- sapply(1:7,function(y) paste0(y + 9,y + 10)) ### for reading-in NCHS data
### Load data 
setwd(paste0(bfolder,'BEdata'))

fname <- paste0('AgeseasoncodaList.RData')
load(fname)

USpopdata <- read.csv("USpop.csv")

agcat <- 5
agernge <- agecatlist[[agcat]][[1]]
agernge[2] <- min(agernge[2],85) 

FSNpopvec <- USpopvec <- NULL
seasoncodaList <- AgeseasoncodaList[[agcat]]

for (k in 1:7) {
  fname <- paste0("NCHS ",seasvec[k] ," population estimates.xls")
  dataset <- read_excel(fname)
  
  stateredvec <- seasoncodaList[[k]][["states"]]
  FSNpop <- sum(dataset[(agernge[1]:agernge[2]) + 1,stateredvec]) ## sums over states/age ranges of relevance 
  FSNpopvec <- c(FSNpopvec,FSNpop)  
  
  year <- eval(parse(text = paste0("20",substr(seasvec[k],1,2))))
  USpop <- sum(USpopdata$pop[which(USpopdata$year==year & (USpopdata$age>=agernge[1] & USpopdata$age<=agernge[2]))])
  USpopvec <- c(USpopvec,USpop)  
}

### example for season 4
seas <- 4


stateredvec <- c("CA","CO") ## only example-not valid


#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
statevec <- unique(as.vector(FSNtestdata$state))
# agcat <- 5
#########################################################################################
#########################################################################################
AgeseasoncodaList <- list()

for (agcat in 1:5) {
  sensdatasel <- list(sensdata[[1]][[agcat]],sensdata[[2]][[agcat]])
  
  cipcr <- sensdatasel[[1]]/100
  seest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
  pcrsens <- c(cipcr[1],seest)
  #########################################################################################
  cirapid <- sensdatasel[[2]]/100
  #########################################################################################
  selograpidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
  lrapidsens <- c(log(cirapid[1]),selograpidest)
  #########################################################################################
  #########################################################################################
  agecatvec <- agecatlist[[agcat]]
  
  seasoncodaList <- list()
  
  for (seas in 1:nseas) { 
    
    stateredvec <- as.vector(unlist(sapply(statevec, function(st) if (sum(FSNtestdata$freq[which(FSNtestdata$state==st &  FSNtestdata$TestedFlu==1 &
                                                                                                 FSNtestdata$died==1 & FSNtestdata$agecat==agcat &
                                                                                                 FSNtestdata$season==seas)]) > 0) st)))
    
    #########################################################################################
    ###  Analyses by state ##################################################################
    #########################################################################################
    codatotarr <- array(0,dim = c(6000,2))
    
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
      ntotvec <- NULL
      testposarr <- array(0,dim = c(2,3))
      #########################################################################################
      ###  Note: These calcuations for outcome P&I   ##########################################
      #########################################################################################
      nttypearr <- array(0,dim = c(2,4))
      
      testposarr <- array(0,dim = c(2,3))
      for (d in dvec){
        selind2 <- which(FSNtestdatasel$died==d)
        ds <- FSNtestdatasel[selind2,]
        
        nttypearr[d+1,1] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==1)]) # PCR
        nttypearr[d+1,2] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==2)]) # RIDT
        nttypearr[d+1,3] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==3)]) # Other/unknown
        nttypearr[d+1,4] <- sum(ds$freq[which(ds$TestedFlu!=1)]) # not tested
        
        ntotvec[d+1] <- nttypearr[d+1,1] + nttypearr[d+1,2] + nttypearr[d+1,3] + nttypearr[d+1,4]
        
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
      #########################################################################################
      #########################################################################################
      selind5 <- which(oshdatsel$season==seas)
      ds3 <- oshdatsel[selind5,]
      
      posh <- sum(ds3$pi[which(ds3$osh==1)])/(sum(ds3$pi))
      
      # FSNfluhosptot <- FSNfluhosp + FSNfludeath ## Non-fatal flu hosp.s
      #########################################################################################
      #########################################################################################
      #########################################################################################
      
      data <- list('FSNfluhosp'=FSNfluhosp,
                   'posh'=posh,'FSNfludeath'=FSNfludeath,
                   'nttypearr'=nttypearr,'testposarr'=testposarr,'pcrsens'=pcrsens,
                   'lrapidsens'=lrapidsens,'ntotvec'=ntotvec)
      
      ptesttot <- t(sapply(1:2, function(k) sum(nttypearr[k,])))
      
      ptestvec10init <- c(nttypearr[,1]/ptesttot)
      ptestvec20init <- c(nttypearr[,2]/ptesttot)
      ptestvec30init <- c(nttypearr[,3]/ptesttot)
      ptestvec40init <- c(nttypearr[,4]/ptesttot)
      
      sens1vecinit <- c(pcrsens[1],pcrsens[1])
      logsens2vecinit <- c(lrapidsens[1],lrapidsens[1])
      sens3vecinit <- c(.4,.4)
      
      sensarrinit <- array(0,dim = c(2,3))
      
      for (k in 1:2){
        sensarrinit[k,1] <- sens1vecinit[k]
        sensarrinit[k,2] <- exp(logsens2vecinit[k])
        sensarrinit[k,3] <- sens3vecinit[k]
      }
      
      ptestarrinit <- array(0,dim = c(2,4))
      
      for (k in 1:2){
        ptestarrinit[k,1] <- ptestvec10init[k]
        ptestarrinit[k,2] <- ptestvec20init[k]
        ptestarrinit[k,3] <- ptestvec30init[k]
        ptestarrinit[k,4] <- ptestvec40init[k]
      }
      
      ptvecinit <- NULL
      
      ptvecinit[1] <- sum(sensarrinit[1,]*ptestarrinit[1,1:3])
      ptvecinit[2] <- sum(sensarrinit[2,]*ptestarrinit[2,1:3])
      
      fluposarrinit <- round(testposarr/sensarrinit)
      
      for (k in 1:2) {
        for (t in 1:3) {
          if (fluposarrinit[k,t] > nttypearr[k,t]) {
            fluposarrinit[k,t] <- nttypearr[k,t]
          }
        }
      }
      
      
      ## Add one to ensure non-zero denominator
      rfluhospinit <- FSNfluhosp /ptvecinit[1]
      rfludeathinit <- FSNfludeath/(1 - posh)/ptvecinit[2]
      
      pfluvecinit <- rowSums(fluposarrinit)/rowSums(nttypearr)
      
      inits <- function(){
        list(
          rfluhosp = rfluhospinit,
          rfludeath = rfludeathinit,
          ptestvec10 = ptestvec10init,
          ptestvec20 = ptestvec20init,
          ptestvec30 = ptestvec30init,
          ptestvec40 = ptestvec40init,
          fluposarr = fluposarrinit,
          pfluvec = pfluvecinit,
          sens1vec = sens1vecinit,
          logsens2vec = logsens2vecinit,
          sens3vec = sens3vecinit
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
      codatotarr <- codatotarr + codaarr
    }
    seasoncodaList[[seas]] <- list(states=stateredvec,codatotarr)
  }

  AgeseasoncodaList[[agcat]] <- seasoncodaList
    #########################################################################################
  #########################################################################################
}
#########################################################################################
setwd(paste0(bfolder,'BEwriteup'))
fname <- paste0('seasoncodaList',agcat,'.RData')
save(seasoncodaList,file = fname)

#########################################################################################
#########################################################################################
