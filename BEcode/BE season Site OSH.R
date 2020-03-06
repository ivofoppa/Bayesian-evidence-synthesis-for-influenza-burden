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
rcodenm <- "Data management Site w testing data.R"
setwd(paste0(bfolder,'BEcode/Data management'))
source(rcodenm)
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
ag <- 4; seas <- 3; st <- "CT"

codaTotList <- list()

for (ag in 1:5) {
  codaseasList <- list()
  for (seas in 1:nseas) {
    codatotls <- array(0, dim = c(niter*3/5,2))
    statels <- stateList[[seas]]
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
          selind2 <- which(Fulldatasel$died==d)
          ds <- Fulldatasel[selind2,]
          
          nttypearr[d+1,1] <- as.numeric(ds$ntest1) # PCR
          nttypearr[d+1,2] <- as.numeric(ds$ntest2) # RIDT
          nttypearr[d+1,3] <- as.numeric(ds$ntest3) # Other/unknown
          nttypearr[d+1,4] <- as.numeric(ds$ntest0) # not tested
          
          ntot[d+1] <- nttypearr[d+1,1] + nttypearr[d+1,2] + nttypearr[d+1,3] + nttypearr[d+1,4]
          
          testposarr[d+1,1] <- as.numeric(ds$testres1) # PCR
          testposarr[d+1,2] <- as.numeric(ds$testres2) # Other/unknown
          testposarr[d+1,3] <- as.numeric(ds$testres3) # Other/unknown
        }
        
        #########################################################################################
        #########################################################################################
        Npideath <- sum(OSHdatasel$pi[which(OSHdatasel$season==seas )])
        oshpideath <- OSHdatasel$pi[which(OSHdatasel$season==seas & OSHdatasel$osh==1)]
        
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
        
        FSNpop <- as.numeric(Fulldatasel$FSNpop)[1]
        
      data <- list('FSNfluhosp'=FSNfluhosp,'FSNpop'=FSNpop,
                   'Npideath'=Npideath,'oshpideath'=oshpideath,'FSNfludeath'=FSNfludeath,
                   'nttypearr'=nttypearr,'testposarr'=testposarr,'pcrsens'=pcrsens,
                   'lrapidsens'=lrapidsens,'ntot'=ntot)
      
      ptesttot <- sapply(1:2, function(k) sum(nttypearr[k,] + 1))
      
      ptest10init <- (nttypearr[,1] + 1)/ptesttot
      ptest20init <- (nttypearr[,2] + 1)/ptesttot
      ptest30init <- (nttypearr[,3] + 1)/ptesttot
      ptest40init <- (nttypearr[,4] + 1)/ptesttot
      
      sens1init <- c(pcrsens[1],pcrsens[1])
      logsens2init <- c(lrapidsens[1],lrapidsens[1])
      sens3init <-c(.4,.4)
      
      sensarrinit <- array(0,dim = c(2,3))
      
        for (k in 1:2){
          sensarrinit[k,1] <- sens1init[k]
          sensarrinit[k,2] <- exp(logsens2init[k])
          sensarrinit[k,3] <- sens3init[k]
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
      
      rfluhospinit <- FSNfluhosp*3/FSNpop
      rfludeathinit <- FSNfludeath/(1 - poshinit)/FSNpop
      
      inits <- function(){
        list(
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
      
      variables <- c('fludeath','fluhosp')
      # variables <- c('USfluhospls')
      # variables <- c('pt')
      # variables <- c('pt')
      
      nadapt <- 10000
      niter <- 10000
      
      model.file <- 'BE season OSH state season.txt'
      
      setwd(paste0(bfolder,'BEmodels'))
      
      j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
      j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
      
      # summary(j.samples)
      
      codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
      codatotls <- codatotls + codals  
    }
    codaseasList[[seas]] <- codatotls
  }
  codaTotList[[ag]] <- codaseasList
}


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
