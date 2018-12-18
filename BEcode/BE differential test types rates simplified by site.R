library(R2jags)
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'
### This code is for one age group, only for season 2014-15.
### The age group has to be selected; separately for the two data sets!
agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

### Age specific estimates from Millman et al. EID 2015, 21 (9):
cipcrlist <- list(c(95.0,82,98.7),c(95.0,82,98.7),c(94.1,81.1,98.7),c(94.1,81.1,98.7),c(86.1,79.6,92.7))
cirapidlist <- list(c(66.7,61.3,71.7),c(66.7,61.3,71.7),c(53.9,47.8,59.8),c(53.9,47.8,59.8),c(20.1,8.8,41.4))
#########################################################################################
### Reading-in data set and making sure DateHosp is actually data variable  #############
#########################################################################################
### Make sure data set is not attached yet; if it is, detach it ...
if (exists('ag'))
{
  rm(ag)
}
for (i in 1:5){
  attached <- search()
  if ('seldata'%in%attached)
  {detach(seldata)}
}

# agcat <- 4

agecatls <- agecatlist[[agcat]]
#########################################################################################
###  Defining influenza seasons #########################################################
#########################################################################################
fromyr <- 2010
toyr <- 2016
fromwk <- 40
towk <- 39

nseas <- toyr - fromyr
### Define 'global N'

mmwrdat0<-read.table("mmwrweeks.txt",header=T)
mmwrdat<-mmwrdat[which((mmwrdat0$wyear==fromyr&mmwrdat0$week>=fromwk)|(mmwrdat0$wyear>fromyr&mmwrdat0$wyear<toyr)|(mmwrdat0$wyear==toyr&mmwrdat0$week<=towk)),]

DVDun <- mmwrdat$dvdweek

N <- length(mmwrdat$dvdweek)
seasbeg <- which(mmwrdat$week==fromwk)
seasend <- which(mmwrdat$week==towk)
time <- (1:N)/N

seaslist <- list()
for (seas in 1:nseas){
  seaslist[[seas]] <- c(max(seasbeg[seas],3),min(seasend[seas],N))
}
#########################################################################################
###  FluSurv-NET data set ###############################################################
#########################################################################################
### Data managment for detection multiplier
setwd(paste0(bfolder,'BEdata'))
dataset <- data.frame(read.csv('20181214_1415burden_bysite.csv'))

### Changing vars into simple values
for (k in 1:length(dataset[1,])) {
  dataset[,k] <- as.vector(dataset[,k])
}

DateHosp <- as.Date(dataset$DateHosp,"%m/%d/%Y")
DateHosp$dataset <- DateHosp

statels <- unique(dataset$State)
#########################################################################################
### NCHS data set: Weekly mortality by age group, week, diagnostic group and ############
### place of death for proportion deaths outside the hospital                ############
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
dataset2 <- data.frame(read.csv('mort2010_16_week_5agstate.csv'))
dataset2$oshospdth <- sapply(dataset2$placdth, function(plce) ifelse(plce==1,1,0) )

dataset2$DVD <- sapply(dataset2$studywk, function(week) mmwrdat0$dvdweek[which(mmwrdat0$studywk==week)])
### Changing vars into simple values
for (k in 1:length(dataset2[1,])) {
  dataset2[,k] <- as.vector(dataset2[,k])
}
stateselind <- which(dataset2$stateoc %in% statels) ## restrict data to states represented in FluSurv-NET data
dataset2 <- dataset2[stateselind,]
###  summing over seasons ###############################################################
dataset2cum <- NULL
for (seas in 1:nseas){
  for (ag in 1:5){
    for (state in statels){
      for (plce in 0:1){
        selind <- which (dataset2$studywk)
      }
    }
  }
}
#########################################################################################
#########################################################################################
#########################################################################################
selindst <- which(dataset$State==state)
datasetstate <- dataset[selindst,]

for (agcat in 1:5){
  selindag <- which(datasetstate$Age >= agecatls[[1]][1] & datasetstate$Age < agecatls[[1]][2])
  datasetagst <- datasetstate[selindag,]
  ### testtypes of influenza positived, deceased or not
  nttype <- array(0,dim = c(2,4))
  ntot <- c(0,0)
  for (d in c(0,1)){
    selind2 <- which(datasetagst$Died==d)
    ntot[d+1] <- length(selind2)
    ds <- datasetagst[selind2,]
    nttype[d+1,1] <- length(which(ds$TestedFlu==1 & ds$TestType==1)) # PCR
    nttype[d+1,2] <- length(which(ds$TestedFlu==1 & ds$TestType==2)) # RIDT
    nttype[d+1,3] <- length(which(ds$TestedFlu==1 & ds$TestType%in%c(3:9))) # Other/unknown
    nttype[d+1,4] <- length(which(ds$TestedFlu!=1)) # Other/unknown
  }
  #########################################################################################
  ###  nttype according to all hosp/died system ###########################################
  #########################################################################################
  nttype[1,] <- colSums(nttype)
  ntot[1] <- sum(ntot)
  #########################################################################################
  cipcr <- cipcrlist[[agcat]]/100
  seest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
  pcrsens <- c(cipcr[1],seest)
  
  testpos <- array(0,dim = c(2,3))
  
  for (d in c(0,1)){
    selind3 <- which(datasetagst$Died==d & datasetagst$TestResult==1)
    ds <- datasetagst[selind3,]
    testpos[d+1,1] <- length(which(ds$TestType==1)) # PCR
    testpos[d+1,2] <- length(which(ds$TestType==2)) # Other/unknown
    testpos[d+1,3] <- length(which(ds$TestType%in%c(3:9))) # Other/unknown
  }
  #########################################################################################
  ###  testpos according to all hosp/died system ##########################################
  #########################################################################################
  testpos[1,] <- colSums(testpos)
  #########################################################################################
  #########################################################################################
  #########################################################################################
  cirapid <- cirapidlist[[agcat]]/100
  
  selograpidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
  lrapidsens <- c(log(cirapid[1]),selograpidest)
  #########################################################################################
  ###  Reding-in input values #############################################################
  #########################################################################################
  inputdata <- read.csv('Burden inputs All Ages.csv',header = T)
  #########################################################################################
  for (col in 1:length(inputdata[1,])){
    inputdata[,col] <- as.vector(inputdata[,col])
    if (col > 2){
      inputdata[,col] <- as.numeric(inputdata[,col])
    }
  }
  ### NEED TO CONTINUE HERE
  posh <- RCdeathoshprop
  
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
