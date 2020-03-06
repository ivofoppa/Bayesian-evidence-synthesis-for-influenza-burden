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
infname <- 'FluSURV-NET-states.RData'
setwd(paste0(bfolder,'BEdata'))
load(infname)
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
agcat <- 5
nseas <- 6

selind <- which(FSNtestdata$agecat==agcat)
FSNtestdata <- FSNtestdata[selind,]

statels <- unique(FSNtestdata$state)

staterev <- sapply(seq_along(FSNtestdata$state), function(x) which(statels==FSNtestdata$state[x]))
FSNtestdata$state <- staterev
nstate <- max(staterev)
#########################################################################################
###  Generating data uations for outcome P&I   ##########################################
#########################################################################################

dataset <- NULL

for (seas in 1:nseas) {
  for (st in 1:nstate) {
    
    selind0 <- which(FSNtestdata$state==st & FSNtestdata$season==seas & FSNtestdata$TestedFlu==0)
    selind1 <- which(FSNtestdata$state==st & FSNtestdata$season==seas & FSNtestdata$TestedFlu!=0)
    ds0 <- FSNtestdata[selind0,]
    ds1 <- FSNtestdata[selind1,]
    
    freq0 <- sum(ds0$freq) # not tested
    freq1 <- sum(ds1$freq) # tested
    if ((freq0 + freq1) > 0) {
      row <- c(seas,st,freq0,freq1)
      dataset <- rbind(dataset,row,deparse.level = 0)
    }
  }
}

dataset <- data.frame(dataset)
colnames(dataset) <- c('season','state','notest','test')
#########################################################################################
#########################################################################################
#########################################################################################
N <- length(dataset$season)

data <- list('tested'=dataset$test, 'tot'=dataset$test + dataset$notest, 'N' = N, 'seas' = dataset$season, 'nseas' = 6)

pinit <- sapply(1:nseas, function(s) sum(dataset$test[which(dataset$season==s)])/
                  (sum(dataset$test[which(dataset$season==s)]) + sum(dataset$notest[which(dataset$season==s)])))

inits <- function(){
  list(
    p = pinit,
    a = rep(1,nseas),
    b = rep(1,nseas)
  )}

# variables <- c('fludeathls')
variables <- c('p')
# variables <- c('pt')
  # variables <- c('pt')
  
nadapt <- 10000
niter <- 10000

model1.str <- "model {
  for (k in 1:N) {

tested[k] ~ dbin(p[seas[k]],tot[k])
}

for (seas in 1:nseas) {
p[seas] ~ dbeta(1,1) 
}
}
"
model2.str <- "model {
  for (k in 1:N) {

tested[k] ~ dbin(p[seas[k]],tot[k])
}

for (seas in 1:nseas) {
p[seas] ~ dbeta(a[seas],b[seas])

a[seas] ~dgamma(0.001,0.001)
b[seas] ~dgamma(0.001,0.001)
}
}
"
model.spec<-textConnection(model2.str)

j.model <- jags.model(file=model.spec,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
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
