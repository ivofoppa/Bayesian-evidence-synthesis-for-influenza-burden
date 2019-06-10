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
ptname <- "20190320_Detection_probs_for Ivo_v2.csv"
ptdata <- read.csv(ptname, header = T)
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
# agcat <- 5
estList <- list()

for (agcat in 1:5) {
  seasestarr <- array(0,dim = c(0,3))
  
  selind1 <- which(ptdata$agegrp5==agcat)
  ptdatsel <- ptdata[selind1,]
  ptdatsel$seas <- 1:7
  
  selind2 <- which(popdata$agecat==agcat)
  agpopsel <- popdata[selind2,]
  
  selind3 <- which(OSHcumdata$agecat==agcat)
  oshdatsel <- OSHcumdata[selind3,]
  
  selind4 <- which(FSNcumdata$agecat==agcat)
  FSNcumdatasel <- FSNcumdata[selind4,]
  
  agecatls <- agecatlist[[agcat]]
  
  FSNflunodeathvec <- FSNcumdatasel$freq[which(FSNcumdatasel$died==0 & FSNcumdatasel$season <= nseas)]
  FSNfludeathvec <- FSNcumdatasel$freq[which(FSNcumdatasel$died==1 & FSNcumdatasel$season <= nseas)]
  
  FSNfluhospvec <- FSNflunodeathvec + FSNfludeathvec
  
  #########################################################################################
  ###  Note: These calcuations for outcome P&I   ##########################################
  #########################################################################################
  ptvec <- NULL
  for (seas in 1:nseas){
    ptseassel <- which(ptdatsel$seas==seas)
    ptvec <- c(ptvec,ptdatsel$prob_d[ptseassel])
  }
  FSNpopvec <- agpopsel$FSNpop
  USpopvec <- agpopsel$USpop
  #########################################################################################
  #########################################################################################
  #########################################################################################
  
  data <- list('FSNfluhospvec'=FSNfluhospvec,'FSNpopvec'=FSNpopvec, 'USpopvec' = USpopvec,
               'ptvec'=ptvec,'nseas'=nseas)
  
  ## Add one to ensure non-zero denominator
  rfluhospvecinit <- FSNfluhospvec /FSNpopvec*3
  
  USfluhospvecinit <- rfluhospvecinit*USpopvec
  
  inits <- function(){
    list(
      rfluhospvec = rfluhospvecinit
    )}
  
  variables <- c('USfluhospvec')
  
  nadapt <- 10000
  niter <- 10000
  
  model.file <- 'BE season check.txt'
  
  setwd(paste0(bfolder,'BEmodels'))
  
  j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  estrow <- quantile(codals,probs = c(.5,.025,.975))
  seasestarr[seas,] <- estrow
}


codalsdiffsimpl <- data.frame(codals)
#########################################################################################
#########################################################################################
