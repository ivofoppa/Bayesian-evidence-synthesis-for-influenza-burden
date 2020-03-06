library(R2jags)
### This code is for one age group, only for season 2014-15.
### The age group has to be selected; separately for the two data sets!
agecatls <- list(c(65,999),'65+')
### Data managment for detection multiplier
bfolder <- 'C:/Users/VOR1/Dropbox/Misc work/Bayesian evidence synthesis/Bayesian-evidence-synthesis-for-influenza-burden/'
setwd(paste0(bfolder,'BEdata'))
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

Pdetect <- 0.32; PdetectSE <- 0.04;
Pdthdetect <- 0.28; PdthdetectSE <- 0.05; RCdeathoshprop <- .669;

data <- list('FSNfluhosp'=FSNfluhosp,'FSNpop'=FSNpop,'Npop'=USpop,
             'posh'=RCdeathoshprop,'FSNfludeath'=FSNfludeath,
             'Pdthdetect'=Pdthdetect,'PdthdetectSE'=PdthdetectSE,
             'Pdetect'=Pdetect,'PdetectSE'=PdetectSE)


ptinit <- c(Pdetect,Pdthdetect)

rfluinit <- FSNfluhosp/FSNpop*3
rfludeathinit <- FSNfludeath/FSNpop*3

pdeathinit <- FSNfludeath/FSNfluhosp

inits <- function(){
  list(
    rfluhosp = rfluinit,
    pdeath = pdeathinit,
    pt = ptinit
    )}

variables <- c('USfluhosp','USfludeath')

model.file <- 'BE monte carlo rates.txt'
nadapt <- 10000
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
