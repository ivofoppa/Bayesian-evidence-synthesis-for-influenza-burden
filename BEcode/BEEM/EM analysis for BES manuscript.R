library(splines)
library(R2jags)
library(stringr)
###############################################################################
###  Ivo M. Foppa, March 2020-----Check for BES paper #########################
###############################################################################
rm(list = ls())
bfolder <- "C:/Users/.../.../" ## project folder

setwd(paste0(bfolder,'...')) ### set to subfolder with data

filename <- paste0('data National 3 ag 2010_18.dat')
datarr <- data.frame(read.table(file = filename,header = T,as.is = T))
###################################################################################################
niter <- 5000
nadapt <- 1000
###################################################################################################
fromyr <- 2012
toyr <- 2017
fromwk <- 40
towk <- 25
###################################################################################################
###################################################################################################
nseas <- toyr-fromyr
###################################################################################################
###############################################################################
################## 7-season analyses ##########################################
###############################################################################
qual1 <- 'Poisson'
model.file <- paste0('model ',qual1,' ',nseas,' seas.txt')
###################################################################################################
DVDun <- unique(datarr$dvd)

seasbeg <- which(substr(DVDun,5,6)==fromwk)
seasend <- which(substr(DVDun,5,6)==towk )

N <- length(DVDun)
for (k in 1:nseas){
  assign(paste("seas",k,sep=""),c(max(seasbeg[k],3),min(seasend[k],N)))
}
###################################################################################################
###  For the dynamic data definition ##############################################################
###################################################################################################
seasls <- sapply(1:nseas, function(seas) paste0('\"seas',seas,'\"=','seas',seas))
seasls2 <- paste0(seasls,collapse = ',')
###################################################################################################
###################################################################################################
fname <- paste0('codaarr check ',qual1," ",niter,' National allcause ',fromyr,'-',substr(toyr,3,4),'.RData')
###################################################################################################
###################################################################################################
variables <- c(sapply(1:nseas, function(seas) paste0('EM',seas,'tot')))
#######################################################################################
variables <- c(sapply(1:nseas, function(seas) paste0('EM',seas,'tot')))

codaarrag <- list()
codaarrtot <- array(0,dim = c(round(niter*3/5),length(variables)))
for (ag in 1:3){
  #######################################################################################
  ###################################################################################################
  agdata <- data.frame(datarr[which(datarr$agcat==ag),])
  
  mort <- agdata$allcause
  pop <- agdata$pop
  AH1 <- agdata$AH1
  AH3 <- agdata$AH3
  B <- agdata$B
  
  #########################################################################################
  #######################################################################################
  N <- length(mort)
  time <- (1:N)/N
  
  ndf <- 4 * nseas
  nsarr <- ns(time,df = ndf)[,]  ## defining basis matrix
  #######################################################################################
  #######################################################################################
  fact <- 10000 ## Flu indicator multiplier
  data <- eval(parse(text = paste0('list("N"=N,"ndf"=ndf,"ns"=nsarr,"mort"=mort,"pop"=pop,"AH3"=AH3*fact,
                                     "AH1"=AH1*fact,"B"=B*fact,',seasls2,')')))
  #######################################################################################
  mod <- lm(mort ~ ns(time, df = ndf))
  smod <- summary(mod)
  coeffls <- as.numeric(smod$coefficients[,1])/mean(pop)
  nsinit <- coeffls[2:(ndf + 1)]
  b0init <- coeffls[1]
  #######################################################################################
  inits <- function() {
    list(
      b0=b0init,
      
      b10=0,
      b11=0,
      b12=0,
      
      b20=0,
      b21=0,
      b22=0,
      
      b30=0,
      b31=0,
      b32=0,
      
      b=nsinit
    )}
  
  setwd(paste0(bfolder,"...")) #### Set to subfolder containing the model
  
  j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
  
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  
  codaarrag[[ag]] <- codals
  codaarrtot <- codaarrtot + codals
  
  cat(paste0('\nAge group ',ag,', National, ',qual1,', toyr=',toyr,': done\n\n'))
  
  setwd(paste0(bfolder,"BEcode/BEEM"))
  save(codaarrag,codaarrtot,model.file, ndf,toyr,file = fname)
}
###################################################################################################
###  Data table ###################################################################################
###################################################################################################
fromyr <- 2012
toyr <- 2017
fromwk <- 40
towk <- 25

nseas <- toyr-fromyr
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
outfile <- paste0('Nation allcause EM check ',qual1," ",fromyr,'-',substr(toyr,3,4),'.txt')
###################################################################################################
###################################################################################################
setwd(paste0(bfolder,"...")) ### Subfolder for writing the resulting file
# load(fname)

seasls <- c(sapply(1:nseas, function(seas) paste0(fromyr + (seas-1),'/',substr(fromyr + seas,3,4))),'All Years')
yrls <- fromyr:toyr
firstline <-paste0('Season\t0-17\t18-64\t65+\tAll Ages')

write.table('',outfile,append=F,row.names=FALSE,col.names=FALSE, quote=FALSE)

### summed over all ages
elarr <- NULL
for (seas in 1:(nseas)){
  elarr <- rbind(elarr,c(seasls[seas],as.vector(quantile(codaarrtot[,seas],prob=c(.5,.05,0.975)))),deparse.level = 0)
}

write.table(elarr,outfile,append=F,row.names=FALSE,col.names=FALSE, quote=FALSE)
###################################################################################################
###################################################################################################