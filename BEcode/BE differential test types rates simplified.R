library(R2jags)
### This code is for one age group, only for season 2014-15.
### The age group has to be selected; separately for the two data sets!
agecatls <- list(c(65,999),'65+')
### Data managment for detection multiplier
bfolder <- 'C:/Users/VOR1/Dropbox/Misc work/Bayesian evidence synthesis/Bayesian-evidence-synthesis-for-influenza-burden/'
setwd(paste0(bfolder,'BEdata'))



dataset <- data.frame(read.csv('FluSurv-Net for sens etc.csv'))

### Changing vars into simple values
for (k in 1:length(dataset[1,])) {
  dataset[,k] <- as.vector(dataset[,k])
}

# dataset$DateHosp <- as.Date(dataset$DateHosp,"%d%B%Y")
# 
# dataset$DateDeath <- as.Date(dataset$DateDeath,"%d%B%Y")
selind <- which(dataset$Age >= agecatls[[1]][1] & dataset$Age < agecatls[[1]][2])
datasetag <- dataset[selind,]
### testtypes of influenza positived, deceased or not
nttype <- array(0,dim = c(2,4))
ntot <- c(0,0)
for (d in c(0,1)){
  selind2 <- which(datasetag$Died==d)
  ntot[d+1] <- length(selind2)
  ds <- datasetag[selind2,]
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

### test sensitivities in 65+ from Millman et al. 2015, adapted by Melissa Rolfes
cipcr <- c(86.1,79.6,92.7)/100
seest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
pcrsens <- c(cipcr[1],seest)

testpos <- array(0,dim = c(2,3))

for (d in c(0,1)){
  selind3 <- which(datasetag$Died==d & datasetag$TestResult==1)
  ds <- datasetag[selind3,]
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
length(which(datasetag$TestedFlu==1))/length(datasetag$PatientNo)
#PCR
cipcr <- c(86.1,79.6,92.7)/100
sepcrest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
pcrsens <- c(cipcr[1],sepcrest)

cirapid <- c(20.1,8.8,41.4)/100
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

agels <- unique(inputdata$ag)

selind <- which(inputdata$ag==agecatls[[2]] & inputdata$season=="2014-15") 
seldata <- inputdata[selind,]
### Make sure data set is not attached yet; if it is, detach it ...
for (i in 1:3){
  attached <- search()
  if ('seldata'%in%attached)
  {detach(seldata)}
}

attach(seldata)

data <- list('FSNfluhosp'=FSNfluhosp,'FSNpop'=FSNpop,'Npop'=USpop,
             'posh'=RCdeathoshprop,'FSNfludeath'=FSNfludeath,
             'nttype'=nttype,'testpos'=testpos,'pcrsens'=pcrsens,'lrapidsens'=lrapidsens,'ntot'=ntot)

pt1init <- nttype[,1]/ntot
pt20init <- nttype[,2]/ntot/(1 - pt1init)
pt30init <- nttype[,3]/ntot/(1 - pt1init - pt20init*(1 - pt1init))

sens1init <- c(pcrsens[1],pcrsens[1])
logsens2init <- c(lrapidsens[1],lrapidsens[1])
sens3init <- c(.4,.4)

pfluinit <- c(testpos[,1]/sens1init + testpos[,2]/exp(logsens2init) + testpos[,3]/sens3init)/rowSums(nttype[,1:3])

fluposinit <- round(pfluinit*nttype[,1:3])

rfluinit <- FSNfluhosp/FSNpop*3

ptestinit <- nttype/rowSums(nttype)
pdeathinit <- FSNfludeath/FSNfluhosp

inits <- function(){
  list(
    rfluhosp = rfluinit,
    pdeath = pdeathinit,
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

nadapt <- 100000
niter <- 100000

model.file <- 'BE differential test types rates simplified.txt'

setwd(paste0(bfolder,'BEmodels'))

j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 

summary(j.samples)

codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
codalsdiffsimpl <- data.frame(codals)

setwd(paste0(bfolder,'BEwriteup'))
save(codalsdiffsimpl,file = 'codalsdiffsimpl.RData')
#########################################################################################
###  ####################################################################################
#########################################################################################
plot(codalsdiffsimpl$USfluhosp,type = 'l')
