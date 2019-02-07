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
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FluSurvdata,sensdata,agseaspop,oshdat
st <- 'NY'
agcat <- 5

selind1 <- which(FluSurvdata$state==st & FluSurvdata$ag==agcat)
FluSurvdatasel <- FluSurvdata[selind1,]

selind2 <- which(agseaspop$state==st & agseaspop$ag==agcat)
agseaspopsel <- agseaspop[selind2,]

selind3 <- which(oshdat$state==st & oshdat$ag==agcat)
oshdatsel <- oshdat[selind3,]

sensdatasel <- list(sensdata[[1]][agcat],sensdata[[2]][agcat])

for (agcat in 1:5){
  agecatls <- agecatlist[[agcat]]
  
  nttypelist <- list()
  ntotlist <- list()
  testposlist <- list()
  poshlist <- list()
  
  for (seas in 1:5){
    nttype <- array(0,dim = c(2,4))
    ntot <- c(0,0)
    for (d in c(0,1)){
      selind2 <- which(FluSurvdatasel$died==d & FluSurvdatasel$season==seas)
      ntot[d+1] <- length(selind2)
      ds <- FluSurvdatasel[selind2,]
      nttype[d+1,1] <- length(which(ds$TestedFlu==1 & ds$TestType==1)) # PCR
      nttype[d+1,2] <- length(which(ds$TestedFlu==1 & ds$TestType==2)) # RIDT
      nttype[d+1,3] <- length(which(ds$TestedFlu==1 & ds$TestType==3)) # Other/unknown
      nttype[d+1,4] <- length(which(ds$TestedFlu!=1)) # not tested
    }
    #########################################################################################
    ###  nttype according to all hosp/died system ###########################################
    #########################################################################################
    nttype[1,] <- colSums(nttype)
    ntot[1] <- sum(ntot)
    
    nttypelist[[seas]] <- nttype
    ntotlist[[seas]] <- ntot
    
    #########################################################################################
    testpos <- array(0,dim = c(2,3))
    for (d in c(0,1)){
      selind3 <- which(FluSurvdatasel$died==d & FluSurvdatasel$TestResult==1 & FluSurvdatasel$season==seas)
      ds <- FluSurvdatasel[selind3,]
      testpos[d+1,1] <- length(which(ds$TestType==1)) # PCR
      testpos[d+1,2] <- length(which(ds$TestType==2)) # Other/unknown
      testpos[d+1,3] <- length(which(ds$TestType==3)) # Other/unknown
    }
    #########################################################################################
    ###  testpos according to all hosp/died system ##########################################
    #########################################################################################
    testpos[1,] <- colSums(testpos)
    testposlist[[seas]] <- testpos
    #########################################################################################
    #########################################################################################
    selind4 <- which(oshdatsel$season==seas)
    ds <- oshdatsel[selind4,]

    poshlist[[seas]] <- ds$freq[which(ds$osh==1)]/(sum(ds$freq))
    
  }
  #########################################################################################
  #########################################################################################
  #########################################################################################
  cipcr <- cipcrlist[[agcat]]/100
  seest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
  pcrsens <- c(cipcr[1],seest)
  
  #########################################################################################
  #########################################################################################
  #########################################################################################
  cirapid <- cirapidlist[[agcat]]/100
  
  selograpidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
  lrapidsens <- c(log(cirapid[1]),selograpidest)
  
  ### NEED TO CONTINUE HERE
    
  data <- list('FSNfluhosp'=FSNfluhosp-FSNfludeath,'FSNpop'=FSNpop,'Npop'=USpop,
               'posh'=posh,'FSNfludeath'=FSNfludeath,
               'nttype'=nttype,'testpos'=testpos,'pcrsens'=pcrsens,'lrapidsens'=lrapidsens,'ntot'=ntot)
  
  pt1init <- sapply(nttype[,1]/ntot, function(x) ifelse(x > .5,min(x,.95),max(x,.1)))
  pt20init <- sapply(nttype[,2]/ntot/(1 - pt1init), function(x) ifelse(x > .5,min(x,.95),max(x,.1)))
  pt30init <- sapply(nttype[,3]/ntot/(1 - pt1init - pt20init*(1 - pt1init)),
                     function(x) ifelse(x > .5,min(x,.95),max(x,.1)))
  
  sens1init <- c(pcrsens[1],pcrsens[1])
  logsens2init <- c(lrapidsens[1],lrapidsens[1])
  sens3init <- c(.4,.4)
  
  pfluinit <- sapply(c(testpos[,1]/exp(logsens1init) + testpos[,2]/exp(logsens2init) + testpos[,3]/sens3init)/rowSums(nttype[,1:3]),
                     function(x) ifelse(x > .5,min(x,.95),max(x,.1)))
  
  fluposinit <- round(pfluinit*nttype[,1:3])
  
  rfluhospinit <- FSNfluhosp/FSNpop*3
  rfludeathinit <- FSNfludeath/FSNpop/(1-posh)
  
  ptestinit <- nttype/rowSums(nttype)
  
  inits <- function(){
    list(
      rfluhosp = rfluhospinit,
      rfludeath = rfludeathinit,
      ptest1 = pt1init,
      ptest20 = pt20init,
      ptest30 = pt30init,
      pflu = pfluinit,
      flupos = fluposinit,
      sens1 = sens1init,
      logsens2 = logsens2init,
      sens3 = sens3init
    )}
  
  variables <- c('USfluhosp','USfludeath')
  # variables <- c('pt')
  # variables <- c('pt')
  
  nadapt <- 10000
  niter <- 10000
  
  model.file <- 'BE differential test types rates simplified.txt'
  
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
