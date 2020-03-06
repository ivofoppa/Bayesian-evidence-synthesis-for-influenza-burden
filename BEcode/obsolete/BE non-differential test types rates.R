library(R2jags)
bfolder <- 'C:/Users/VOR1/Dropbox/Misc work/Bayesian evidence synthesis/Bayesian-evidence-synthesis-for-influenza-burden/'
### This code is for one age group, only for season 2014-15.
### The age group has to be selected; separately for the two data sets!
agecatlist <- list(list(c(0,5),'<5'),
                   list(c(5,18),'5 to 17'),
                   list(c(18,50),'18-49'),
                   list(c(50,65),'50-64'),
                   list(c(65,999),'65+'))

### Age specific estimates from Millman et al. EID 2015, 21 (9):
cipcrlist <- list(c(95.0,82,98.7),c(95.0,82,98.7),c(94.1,81.1,98.7),c(94.1,81.1,98.7),c(86.1,79.6,92.7))
cirapidlist <- list(c(66.7,61.3,71.7),c(66.7,61.3,71.7),c(53.9,47.8,59.8),c(53.9,47.8,59.8),c(20.1,8.8,41.4))
# agcat <- 2
#########################################################################################
### Function to conditionally remove objects; copied from ###############################
### https://stackoverflow.com/questions/7172568/write-a-function-to-remove-object-if-it-exists
#########################################################################################
ifrm <- function(ls, env = globalenv()) {
  for (el in ls){
    obj <- deparse(substitute(el))
    if(exists(obj, envir = env)) {
      rm(list = obj, envir = env)
    }
  }
}
#########################################################################################
for (agcat in 1:5){
  agecatls <- agecatlist[[agcat]]
  
  ### Data managment for detection multiplier
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
cipcr <- cipcrlist[[agcat]]/100
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

agels <- unique(inputdata$ag)

selind <- which(inputdata$ag==agecatls[[2]] & inputdata$season=="2014-15") 
seldata <- inputdata[selind,]
#########################################################################################
### Make sure data set is not attached yet; if it is, detach it ...

ifrm(c(ag,Pdetect,PdetectSE,Pdthdetect,PdthdetectSE,RCdeathoshprop))

for (i in 1:5){
  attached <- search()
  if ('seldata'%in%attached)
  {detach(seldata)}
}
#########################################################################################
attach(seldata)

posh <- RCdeathoshprop

data <- list('FSNfluhosp'=FSNfluhosp,'FSNpop'=FSNpop,'Npop'=USpop,
             'posh'=posh,'FSNfludeath'=FSNfludeath,
             'nttype'=nttype,'testpos'=testpos,'pcrsens'=pcrsens,'lrapidsens'=lrapidsens,'ntot'=ntot)

pt10init <- sapply(nttype[,1]/ntot, function(x) ifelse(x > .5,min(x,.95),max(x,.1)))
pt20init <- sapply(nttype[,2]/ntot/(1 - pt10init), function(x) ifelse(x > .5,min(x,.95),max(x,.1)))
pt30init <- sapply(nttype[,3]/ntot/(1 - pt10init - pt20init*(1 - pt10init)),
                   function(x) ifelse(x > .5,min(x,.95),max(x,.1)))

sens1init <- pcrsens[1]
logsens2init <- lrapidsens[1]
sens3init <- c(.4,.4)

pfluinit <- c(testpos[,1]/sapply(nttype[,1], function(x) max(x,1))/sens1init,testpos[,2]/sapply(nttype[,2], function(x) max(x,1))/exp(logsens2init),testpos[,3]/sapply(nttype[,3], function(x) max(x,1))/sens3init,(testpos[,1]/sens1init + testpos[,2]/exp(logsens2init) + testpos[,3]/sens3init)/sapply(rowSums(nttype[,1:3]),function(x) max(x,1))
                  ,deparse.level = 1)+.01

pfluinit <- array(sapply(pfluinit, function(x) ifelse(x==0,.2,ifelse(x>1,.9,x))),dim = c(2,4))

fluposinit <- round(pfluinit[,1:3]*nttype[,1:3])


inits <- function(){
  list(
    rfluhosp = rfluhospinit,
    rfludeath = rfludeathinit,
    ptest10 = pt1init,
    ptest20 = pt20init,
    ptest30 = pt30init,
    pflu = pfluinit,
    flupos = fluposinit,
    sens1 = sens1init,
    logsens2 = logsens2init,
    sens3 = sens3init
  )}


variables <- c('USfluhosp','USfludeath')

model.file <- 'BE non-differential test types rates.txt'
niter <- 10000
nadapt <- 10000
setwd(paste0(bfolder,'BEmodels'))

j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 

summary(j.samples)

codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
codals <- data.frame(codals)

ls <- codals[,1]
plot(ls,ylim = c(0,max(ls)),type = 'l')

