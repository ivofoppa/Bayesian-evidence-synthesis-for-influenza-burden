rm(list = ls())

library(R2jags)
#########################################################################################
#########################################################################################
bfolder <- 'C:/Users/VOR1/Dropbox/Misc work/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

### Load data 
# infname <- 'FluSURV-NET-site-osh-7seas.RData'
# rcodenm <- "Data management Site w testing data.R"
# setwd(paste0(bfolder,'BEcode/Data management'))
# source(rcodenm)
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
### Load data 
infname <- 'FluSURV-NET-site-osh-7seas.RData'
load(infname)
#########################################################################################
# ag <- 2; seas <- 3; statels <- stateList[[seas]]; st <- statels[3]
nseas <- 7
testTotList <- list()

for (ag in 1:5) {
  testseasList <- list()
  for (seas in 1:nseas) {
    statels <- stateList[[seas]]
    TestList <- list()
    for (st in statels) {
      selind1 <- which(Fulldata$agcat==ag & Fulldata$state==st & Fulldata$season==seas)
      Fulldatasel <- Fulldata[selind1,]
      
      selind2 <- which(OSHdata$agcat==ag & OSHdata$state==st & OSHdata$season==seas)
      OSHdatasel <- OSHdata[selind2,]
      
      sensdatasel <- list(sensdata[[1]][[ag]],sensdata[[2]][[ag]])
      
      agecatls <- agecatlist[[ag]]
      
      FSNfluhosp <- as.numeric(Fulldatasel$Noutcome[which(Fulldatasel$died==0)])
      FSNfludeath <- as.numeric(Fulldatasel$Noutcome[which(Fulldatasel$died==1 & Fulldatasel$season <= nseas)])
      
      Npideath <- oshpideath <- NULL
      
      #########################################################################################
      ###  Note: These calculations for outcome P&I   #########################################
      #########################################################################################
      nttypearr <- array(0,dim = c(2,4))
      testposarr <- array(0,dim = c(2,3))
      ntot <- NULL
      
      for (d in c(0,1)){
        selind3 <- which(Fulldatasel$died==d)
        if (length(selind3) > 0) {
          ds <- Fulldatasel[selind3,]
          
          nttypearr[d+1,1] <- as.numeric(ds$ntest1) # PCR
          nttypearr[d+1,2] <- as.numeric(ds$ntest2) # RIDT
          nttypearr[d+1,3] <- as.numeric(ds$ntest3) # Other/unknown
          nttypearr[d+1,4] <- as.numeric(ds$ntest0) # not tested
          
          ntot[d+1] <- nttypearr[d+1,1] + nttypearr[d+1,2] + nttypearr[d+1,3] + nttypearr[d+1,4]
          
          testposarr[d+1,1] <- as.numeric(ds$testres1) # PCR
          testposarr[d+1,2] <- as.numeric(ds$testres2) # Other/unknown
          testposarr[d+1,3] <- as.numeric(ds$testres3) # Other/unknown
          
        } else {
          nttypearr[d+1,] <- c(0,0,0,0)
          testposarr[d+1,] <- c(0,0,0)
        }
      }
      nottestls <- sapply(nttypearr[,4]/sapply(rowSums(nttypearr),function(x) max(1,x)),function(t) ifelse(t==0,NA,1 - t))
      testposls <- testposarr/nttypearr[,1:3]
      tlist <- list(nottestls,testposls)
      TestList[[st]] <- tlist
    }
    testseasList[[seas]] <- TestList
  }
  testTotList[[ag]] <- testseasList
}
#########################################################################################
###  Computing the number of tests performed (all tests)   ##############################
#########################################################################################
nseas <- 7
numtestTotList <- list()

for (ag in 1:5) {
  testseasList <- list()
  for (seas in 1:nseas) {
    statels <- stateList[[seas]]
    TestList <- list()
    for (st in statels) {
      selind1 <- which(Fulldata$agcat==ag & Fulldata$state==st & Fulldata$season==seas)
      Fulldatasel <- Fulldata[selind1,]
      
      selind2 <- which(OSHdata$agcat==ag & OSHdata$state==st & OSHdata$season==seas)
      OSHdatasel <- OSHdata[selind2,]
      
      sensdatasel <- list(sensdata[[1]][[ag]],sensdata[[2]][[ag]])
      
      agecatls <- agecatlist[[ag]]
      
      FSNfluhosp <- as.numeric(Fulldatasel$Noutcome[which(Fulldatasel$died==0)])
      FSNfludeath <- as.numeric(Fulldatasel$Noutcome[which(Fulldatasel$died==1 & Fulldatasel$season <= nseas)])
      
      Npideath <- oshpideath <- NULL
      
      #########################################################################################
      ###  Note: These calculations for outcome P&I   #########################################
      #########################################################################################
      nttypearr <- array(0,dim = c(2,4))
      testposarr <- array(0,dim = c(2,3))
      ntot <- NULL
      
      for (d in c(0,1)){
        selind3 <- which(Fulldatasel$died==d)
        if (length(selind3) > 0) {
          ds <- Fulldatasel[selind3,]
          
          nttypearr[d+1,1] <- as.numeric(ds$ntest1) # PCR
          nttypearr[d+1,2] <- as.numeric(ds$ntest2) # RIDT
          nttypearr[d+1,3] <- as.numeric(ds$ntest3) # Other/unknown
          nttypearr[d+1,4] <- as.numeric(ds$ntest0) # not tested
          
          ntot[d+1] <- nttypearr[d+1,1] + nttypearr[d+1,2] + nttypearr[d+1,3] 
          
          testposarr[d+1,1] <- as.numeric(ds$testres1) # PCR
          testposarr[d+1,2] <- as.numeric(ds$testres2) # Other/unknown
          testposarr[d+1,3] <- as.numeric(ds$testres3) # Other/unknown
          
        } else {
          nttypearr[d+1,] <- c(0,0,0,0)
          testposarr[d+1,] <- c(0,0,0)
          ntot[d+1] <- 0
        }
      }
      TestList[[st]] <- ntot
    }
    testseasList[[seas]] <- TestList
  }
  numtestTotList[[ag]] <- testseasList
}
