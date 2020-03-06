library(R2jags)
#########################################################################################
#########################################################################################
bfolder <- 'C:/Users/VOR1/Dropbox/Misc work/Bayesian-evidence-synthesis-for-influenza-burden/'
### Load data 
# infname <- 'FluSURV-NET-site-osh-7seas.RData'
rcodenm <- "Data management Nation w testing data.R"
setwd(paste0(bfolder,'BEcode/Data management'))
source(rcodenm)

rm(list = ls())

bfolder <- 'C:/Users/VOR1/Dropbox/Misc work/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

setwd(paste0(bfolder,'BEdata'))
infname <- 'FluSURV-NET-Nation-osh-7seas.RData'
load(infname)
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
# ag <- 5; 
nadapt <- 10000
niter <- 10000

model.file <- 'BE season OSH beta dirichlet lognorm sens test.txt'

Fulldata <- Fulldata[which(Fulldata$season >= 3),]
Fulldata$season <- Fulldata$season - 2

OSHcumdata <- OSHcumdata[which(OSHcumdata$season >= 3),]
OSHcumdata$season <- OSHcumdata$season - 2
nseas <- max(OSHcumdata$season)

codaagList <- list()
# codatotls <- array(0,dim = c(niter*3/5, nseas*2))

# spec <- "ptarr"
spec <- "sensarr"

outfile <- paste0("codafileFSN_2012-17_",spec,".RData")
seasageprop <- array(0,dim = c(nseas*5,2))

for (seas in 1:nseas) {
  selind <- which(Fulldata$season==seas)
  datapart <- Fulldata[selind,]
  popsum <- 0
  popls <- NULL
  for (ag in 1:5) {
    selind2 <- which(datapart$agcat==ag)
    pop <- datapart$FSNpop[selind2][1]
    popsum <- popsum + pop
    popls <- c(popls,pop)
  }
  popls <- popls/popsum
  seasageprop[((1:5)+(seas-1)*5),1] <- 1:5
  seasageprop[((1:5)+(seas-1)*5),2] <- popls
}

for (ag in 1:5) {
  selind1 <- which(Fulldata$agcat==ag )
  Fulldatasel <- Fulldata[selind1,]
  
  selind2 <- which(OSHcumdata$agcat==ag)
  OSHcumdatasel <- OSHcumdata[selind2,]
  
  sensdatasel <- list(sensdata[[1]][[ag]],sensdata[[2]][[ag]])
  
  agecatls <- agecatlist[[ag]]
  
  FSNfluhospnonfatls <- sapply(1:nseas, function(s) {
    selind <- which(Fulldatasel$season==s & Fulldatasel$died==0)
    ifelse(length(selind)>0,sum(as.numeric(Fulldatasel$Noutcome[selind])))
  })
  
  FSNfluhospfatls <-  sapply(1:nseas, function(s) {
    selind <- which(Fulldatasel$died==1 & 
                      Fulldatasel$season==s)
    
    ifelse(length(selind)>0,sum(as.numeric(Fulldatasel$Noutcome[selind])),0)
  })
  
  Npideathls <- oshpideathls <- NULL
  nttypearr <- array(0,dim = c(2,4,nseas))
  ntotarr <- array(0,dim = c(2,nseas))
  testposarr <- array(0,dim = c(2,3,nseas))
  poshls <- NULL
  #########################################################################################
  ###  Note: These calculations for outcome P&I   #########################################
  #########################################################################################
  for (seas in 1:nseas){
    nttype <- array(0,dim = c(2,4))
    ntot <- c(0,0)
    testpos <- array(0,dim = c(2,3))
    for (d in c(0,1)){
      selind2 <- which(Fulldatasel$died==d & Fulldatasel$season==seas)
      if (length(selind2) > 0) {
        ds <- Fulldatasel[selind2,]
        
        nttype[d+1,1] <- ds$ntest1 # PCR
        nttype[d+1,2] <- ds$ntest2 # RIDT
        nttype[d+1,3] <- ds$ntest3 # Other/unknown
        nttype[d+1,4] <- ds$ntest0 # not tested
        
        ntot[d+1] <- nttype[d+1,1] + nttype[d+1,2] + nttype[d+1,3] + nttype[d+1,4]
        
        testpos[d+1,1] <- ds$testres1 # PCR
        testpos[d+1,2] <- ds$testres2 # Other/unknown
        testpos[d+1,3] <- ds$testres3 # Other/unknown
      } else {
        nttype[d+1,1] <- 0 # PCR
        nttype[d+1,2] <- 0 # RIDT
        nttype[d+1,3] <- 0 # Other/unknown
        nttype[d+1,4] <- 0 # not tested
        
        ntot[d+1] <- 0
        
        testpos[d+1,1] <- 0 # PCR
        testpos[d+1,2] <- 0 # Other/unknown
        testpos[d+1,3] <- 0 # Other/unknown
      }
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
    Npideathls[seas] <- sum(OSHcumdatasel$pi[which(OSHcumdatasel$season==seas )])
    oshpideathls[seas] <- OSHcumdatasel$pi[which(OSHcumdatasel$season==seas & OSHcumdatasel$osh==1)]
  }
  
  #########################################################################################
  #########################################################################################
  cipcr <- sensdatasel[[1]]/100
  lsepcrest <- ((log(cipcr[3]) - log(cipcr[1])) + (log(cipcr[1]) - log(cipcr[2])))/2/1.96
  lpcrsens <- c(log(cipcr[1]),lsepcrest)
  #########################################################################################
  #########################################################################################
  #########################################################################################
  cirapid <- sensdatasel[[2]]/100
  lserapidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
  lrapidsens <- c(log(cirapid[1]),lserapidest)
  
  FSNpopls <- sapply(1:nseas, function(s) Fulldatasel$FSNpop[which(Fulldatasel$season==s)][1])
  USpopls  <- sapply(1:nseas, function(s) USpopdata$USpop[which(USpopdata$agcat==ag & USpopdata$season==s)])
  
  alpha <- c(1,1,1,1) ### Dirichlet prior
  
  data <- list('FSNfluhospnonfatls'=FSNfluhospnonfatls,'FSNfluhospfatls'=FSNfluhospfatls,
               'FSNpopls'=FSNpopls, 'USpopls' = USpopls,
               'Npideathls'=Npideathls,'oshpideathls'=oshpideathls,
               'nttypearr'=nttypearr,'testposarr'=testposarr,'lpcrsens'=lpcrsens,
               'lrapidsens'=lrapidsens,'ntotarr'=ntotarr,'alpha'=alpha,'nseas'=nseas)
  
  ptesttot <- t(sapply(1:2, function(k) sapply(1:nseas, function(s) sum(nttypearr[k,,s] +1))))
  
  FSNfluhospls <- FSNfluhospfatls + FSNfluhospnonfatls
  
  pdinit <- FSNfluhospfatls/FSNfluhospls ## initial value of prob. hosp. death
  
  ptestarrinit <- array(0,dim = c(2,4,nseas))
  
  for (t in 1:4) {
    ptestarrinit[,t,] <- (nttypearr[,t,1:nseas] + 1)/ptesttot
  }
  
  sens1arrinit <- sapply(1:nseas, function(s) exp(c(lpcrsens[1],lpcrsens[1])))
  sens2arrinit <- sapply(1:nseas, function(s) exp(c(lrapidsens[1],lrapidsens[1])))
  sens3arrinit <- sapply(1:nseas, function(s) c(.4,.4))
  
  sensarrinit <- array(0,dim = c(2,3,nseas))
  
  for (seas in 1:nseas) {
    for (k in 1:2){
      sensarrinit[k,1,seas] <- sens1arrinit[k,seas]
      sensarrinit[k,2,seas] <- sens2arrinit[k,seas]
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
  
  ptarrinit <- ptestarrinit[,1,]*sensarrinit[,1,] + ptestarrinit[,2,]*sensarrinit[,2,] + 
    ptestarrinit[,3,]*sensarrinit[,3,]
  
  ## Add one to ensure non-zero denominator
  poshlsinit <- oshpideathls/Npideathls
  
  rfluhosplsinit <- (FSNfluhospnonfatls/ptarrinit[1,] + FSNfluhospfatls /ptarrinit[2,]) /FSNpopls
  
  fludeathlsinit <- round(rfluhosplsinit*FSNpopls*pdinit/(1 - poshlsinit))
  
  rfludeathlsinit <- FSNfluhospfatls/FSNpopls/(1 - poshlsinit)
  fluhosplsinit <- round(rfluhosplsinit*FSNpopls)
  
  fluhospfatlsinit <- round(FSNfluhospfatls/ptarrinit[2,])
  
  USfludeathlsinit <- round(rfludeathlsinit*USpopls)
  USfluhosplsinit <- round(rfluhosplsinit*USpopls)
  
  inits <- function(){
    list(
      rfluhospls = rfluhosplsinit,
      fluhospls = fluhosplsinit,
      fluhospfatls = fluhospfatlsinit,
      fludeathls = fludeathlsinit,
      ptestarr = ptestarrinit,
      fluposarr = fluposarrinit,
      sens1arr = sens1arrinit,
      sens2arr = sens2arrinit,
      sens3arr = sens3arrinit,
      poshls = poshlsinit,
      USfludeathls = USfludeathlsinit,
      USfluhospls = USfluhosplsinit
    )}
  
  dim <- length(eval(parse(text = paste0("dim(",spec,"init)"))))
  # variables <- c(paste0(spec,"[1:2,",ifelse(dim==3,"1:3,",""),"1:5]"))
  # variables <- c(paste0("poshls[1:",nseas,"]"))
  # variables <- c(paste0("pd[1:",nseas,"]"))
  # variables <- c(paste0("ptestarr[1:2,1:4,1:",nseas,"]"))
  
  setwd(paste0(bfolder,'BEmodels'))
  
  j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
  
  cat("\nAge group ",ag," done!!!\n\n")
  
  # summary(j.samples)
  
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  
  codaagList[[ag]] <- codals
  # codatotls <- codatotls + codals  
}
setwd(paste0(bfolder,'BEmcmc'))

save(codaagList,sensdata,file = outfile)
###################################################################################################
###################################################################################################
codals <- codaagList[[5]]
poshls <- sapply(seq(1,nseas,1),function(seas) median(codals[,seas]))
ptarrls <- array(sapply(seq(1,nseas*2,1),function(seas) median(codals[,seas])),dim = c(2,5))

sensarrls <- array(sapply(seq(1,nseas*2*3,1),function(seas) median(codals[,seas])),dim = c(2,3,5))

ptestarrls <- array(sapply(seq(1,nseas*2*4,1),function(seas) median(codals[,seas])),dim = c(2,4,5))

variables
