library(R2jags)
### Data managment for detection multiplier
bfolder <- 'C:/Users/VOR1/Dropbox/Misc work/Bayesian evidence synthesis/BEanalysis/BE Git project/'
setwd(paste0(bfolder,'BEdata'))

dataset <- data.frame(read.csv('FluSurv-Net for sens etc.csv'))

### Changing vars into simple values
for (k in 1:length(dataset[1,])) {
  dataset[,k] <- as.vector(dataset[,k])
}

# dataset$DateHosp <- as.Date(dataset$DateHosp,"%d%B%Y")
# 
# dataset$DateDeath <- as.Date(dataset$DateDeath,"%d%B%Y")
delind <- which(dataset$Age >= 65)
dataset65 <- dataset[delind,]
### testtypes of influenza positived, deceased or not
nttype <- array(0,dim = c(2,4))
ntot <- c(0,0)
for (d in c(0,1)){
  selind2 <- which(dataset65$Died==d)
  ntot[d+1] <- length(selind2)
  ds <- dataset65[selind2,]
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
  selind3 <- which(dataset65$Died==d & dataset65$TestResult==1)
  ds <- dataset65[selind3,]
  testpos[d+1,1] <- length(which(ds$TestType==1)) # PCR
  testpos[d+1,2] <- length(which(ds$TestType==2)) # Other/unknown
  testpos[d+1,3] <- length(which(ds$TestType%in%c(3:9))) # Other/unknown
}
#########################################################################################
###  testpos according to all hosp/died system ##########################################
#########################################################################################
testpos[1,] <- colSums(testpos)
#########################################################################################
#PCR
cipcr <- c(86.1,79.6,92.7)/100
sepcrest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
pcrsens <- c(cipcr[1],sepcrest)

cirapid <- c(20.1,8.8,41.4)/100
selograpidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
lograpidsens <- c(log(cirapid[1]),selograpidest)

#########################################################################################
###  Reding-in other input values #######################################################
#########################################################################################
Assign <- function(Names, Values) {
  sapply(seq_along(Names), function(i){assign(Names[i], Values[i], envir=.GlobalEnv)});
  invisible()
}
#########################################################################################

varnames <- c('Npop','FSNpop','FSNfluhosp','fludeath','RCdeath','RCdeathosh','perctested',
              'perctestedSE','sens','sensSE','Pdetect','PdetectSE','Multiplier','perctesteddth',
              'perctesteddthSE','sensdth','sensdthSE','C-HRatio','MAprop','MAprop-Min','MAprop-Max',
              'Pdthdetect','PdthdetectSE')

varvalues <- c(46243211,3521076,10864,405,880220,588976,0.51903,0.04131,0.61674,0.05205,0.32011,
               0.03713,3.12395,0.4323,0.0559,0.64,0.09337,11,0.56,0.53,0.6,0.28,0.05)


Assign(varnames,varvalues)

nfluhospobs <- c(FSNfluhosp,fludeath)

RCdeathish <- RCdeath-RCdeathosh

Pdetect <- c(Pdetect,Pdthdetect)
PdetectSE <- c(PdetectSE,PdthdetectSE) 

data <- list('nfluhospobs'=nfluhospobs,'FSNpop'=FSNpop,'Npop'=Npop,
             'RCdeath'=RCdeath,'RCdeathish'=RCdeathish,
             'Pdetect'=Pdetect,'PdetectSE'=PdetectSE)

rfluinit <- nfluhospobs[1]/FSNpop*3
pdihinit <- RCdeathish/RCdeath
phospinit <- .9

pdeathinit <- 0.1

inits <- function(){
  list(
    rflu = rfluinit,
    phosp = phospinit,
    pdeath = pdeathinit,
    pdih = pdihinit
  )}

variables <- c('USfluhosp','USfludeath')

model.file <- 'BE monte carlo rates.txt'
nadapt <- 100000
niter <- 100000
setwd(paste0(bfolder,'BEmodels'))

j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 

summary(j.samples)

codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
codalsMC <- data.frame(codals)

setwd(paste0(bfolder,'BEwriteup'))
save(codalsMC,file = 'codalsMC.RData')

# ls <- codals[,1]
# plot(ls,ylim = c(0,max(ls)),type = 'l')
# tail(codals)
