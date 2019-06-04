library(R2jags)
#########################################################################################
#########################################################################################
rm(list = ls())
# agcat <- 5
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
statevec <- unique(as.vector(FSNtestdata$state))

AgeseasoncodaList <- list()

for (agcat in 1:5) {
  
  cat("\nRunning age cat ",agcat,"...\n\n")
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
    
    stateredvec <- as.vector(unlist(sapply(statevec, function(st) if 
                                           (sum(FSNtestdata$freq[which(FSNtestdata$state==st &  FSNtestdata$TestedFlu==1 & 
                                                                       FSNtestdata$agecat==agcat & 
                                                                       FSNtestdata$season==seas)]) > 0) st)))
    
    #########################################################################################
    ###  Analyses by state ##################################################################
    #########################################################################################
    codatotarr <- array(0,dim = c(6000,1))
    
    for (st in stateredvec) {
      selind1 <- which(FSNtestdata$agecat==agcat & FSNtestdata$state==st & FSNtestdata$season==seas)
      FSNtestdatasel <- data.frame(FSNtestdata[selind1,])
      
      selind4 <- which(FSNdata$agecat==agcat & FSNdata$state==st & FSNdata$season==seas)
      FSNdatasel <- FSNdata[selind4,]
      
      FSNfluhosp <- sum(FSNdatasel$freq)
      
      nttypevec <- NULL
      ntotvec <- NULL
      testposvec <- NULL
      #########################################################################################
      ###  Note: These calcuations for outcome P&I   ##########################################
      #########################################################################################

      testposvec <- NULL
      ds <- FSNtestdatasel
      
      nttypevec[1] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==1)]) # PCR
      nttypevec[2] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==2)]) # RIDT
      nttypevec[3] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==3)]) # Other/unknown
      nttypevec[4] <- sum(ds$freq[which(ds$TestedFlu!=1)]) # not tested
      
      ntot <- nttypevec[1] + nttypevec[2] + nttypevec[3] + nttypevec[4]
      
      testposvec[1] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==1 & ds$TestResult==1)]) # PCR
      testposvec[2] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==2 & ds$TestResult==1)]) # Other/unknown
      testposvec[3] <- sum(ds$freq[which(ds$TestedFlu==1 & ds$TestType==3 & ds$TestResult==1)]) # Other/unknown
      
      #########################################################################################
      ###  testpos according to all hosp/died system ##########################################
      #########################################################################################
      ### Continue here: the analysis is by season and state, so all data has to be prepared only for one season
      #########################################################################################

      data <- list('FSNfluhosp'=FSNfluhosp,
                   'nttypevec'=nttypevec,'testposvec'=testposvec,'pcrsens'=pcrsens,
                   'lrapidsens'=lrapidsens,'ntot'=ntot)
      
      ptesttot <- sum(nttypevec)
      
      ptest10init <- c(nttypevec[1]/ptesttot)
      ptest20init <- c(nttypevec[2]/ptesttot)
      ptest30init <- c(nttypevec[3]/ptesttot)
      ptest40init <- c(nttypevec[4]/ptesttot)
      
      sens1init <- pcrsens[1]
      logsens2init <- lrapidsens[1]
      sens3init <- .4
      
      sensvecinit <- NULL
      
        sensvecinit[1] <- sens1init
        sensvecinit[2] <- exp(logsens2init)
        sensvecinit[3] <- sens3init
      
      ptestvecinit <- NULL
      
        ptestvecinit[1] <- ptest10init
        ptestvecinit[2] <- ptest20init
        ptestvecinit[3] <- ptest30init
        ptestvecinit[4] <- ptest40init

      ptinit <- sum(sensvecinit*ptestvecinit[1:3])

      fluposvecinit <- round(testposvec/sensvecinit)
      
         for (t in 1:3) {
          if (fluposvecinit[t] > nttypevec[t]) {
            fluposvecinit[t] <- nttypevec[t]
          }
      }
      
      
      ## Add one to ensure non-zero denominator
      rfluhospinit <- FSNfluhosp /ptinit

      pfluinit <- sum(fluposvecinit)/sum(nttypevec)
      
      inits <- function(){
        list(
          rfluhosp = rfluhospinit,
           ptest10 = ptest10init,
          ptest20 = ptest20init,
          ptest30 = ptest30init,
          ptest40 = ptest40init,
          fluposvec = fluposvecinit,
          pflu = pfluinit,
          sens1 = sens1init,
          logsens2 = logsens2init,
          sens3 = sens3init
        )}
      
      # variables <- c('fludeath')
      variables <- c('fluhosp')
      # variables <- c('pt')
      # variables <- c('pt')
      
      nadapt <- 10000
      niter <- 10000
      
      model.file <- 'BE season state hosp.txt'
      
      setwd(paste0(bfolder,'BEmodels'))
      
      j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
      j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
      
      codaarr <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
      codatotarr <- codatotarr + codaarr
    }
    cat("Season ",seas)
    
    seasoncodaList[[seas]] <- list(states=stateredvec,codatotarr)
  }
  assign(paste0("seasoncodaList",agcat),seasoncodaList)
    #########################################################################################
  #########################################################################################
}
#########################################################################################
AgeseasoncodaList <- list(seasoncodaList1,seasoncodaList2,seasoncodaList3,seasoncodaList4,seasoncodaList5)

setwd(paste0(bfolder,'BEdata'))
fname <- paste0('AgeseasoncodaList_Hosp.RData')
save(AgeseasoncodaList,nseas,file = fname)

#########################################################################################
#########################################################################################
