library(R2jags)

Ntest <- 100
Ntestpos <- 99

pcrsens <- c(0.7,0.03)


Nfluposinit <- round(Ntestpos/pcrsens[1])
sensinit <- pcrsens[1]

Ntestposls <- 1:100
sensls <- NULL

for (Ntestpos in Ntestposls) {
  data <- list('Ntestpos'=Ntestpos,'Ntest'=Ntest, 'pcrsens' = pcrsens)
  
  inits <- function(){
    list(
      Nflupos = Nfluposinit,
      sens = sensinit
    )}
  
  variables <- c('sens')
  
  nadapt <- 10000
  niter <- 10000
  model.file <- 'Sensitivity model.txt'
  
  setwd(paste0(bfolder,'BEmodels'))
  
  j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3,quiet = T)
  j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5,quiet = T) 
  
  sm <- summary(j.samples)
  sensls <- c(sensls,as.numeric(sm$statistics["Mean"]))
}

plot(Ntestposls,sensls,type = 'l',ylim = c(0,1),ylab ="Posterior mean of sens.", xlab = "# of pos. tests")
lines(Ntestposls,rep(pcrsens[1],length(Ntestposls)),col = 'red')
lines(Ntestposls,sensls)
