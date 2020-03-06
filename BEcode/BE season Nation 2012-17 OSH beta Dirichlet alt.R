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
nseas <- 7
#########################################################################################
#########################################################################################
###  FluSurv-NET data set ###############################################################
###  Subsetting by site and age group ###################################################
#########################################################################################
###  FSNdata,sensdata,agseaspop,oshdat
# ag <- 5; 
nadapt <- 10000
niter <- 10000

model.file <- 'BE season OSH beta dirichlet alt.txt'

Fulldata <- Fulldata[which(Fulldata$season >= 3),]
Fulldata$season <- Fulldata$season - 2

OSHcumdata <- OSHcumdata[which(OSHcumdata$season >= 3),]
OSHcumdata$season <- OSHcumdata$season - 2
nseas <- max(OSHcumdata$season)

codaagList <- list()
codatotls <- array(0,dim = c(niter*3/5, nseas*2))

outfile <- "codafileFSN_OSH_beta_dirichlet_alt_2012-17.RData"

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
  seest <- ((cipcr[3] - cipcr[1]) + (cipcr[1] - cipcr[2]))/2/1.96
  pcrsens <- c(cipcr[1],seest)
  #########################################################################################
  #########################################################################################
  #########################################################################################
  cirapid <- sensdatasel[[2]]/100
  
  selograpidest <- ((log(cirapid[3]) - log(cirapid[1])) + (log(cirapid[1]) - log(cirapid[2])))/2/1.96
  lrapidsens <- c(log(cirapid[1]),selograpidest)
  
  FSNpopls <- sapply(1:nseas, function(s) Fulldatasel$FSNpop[which(Fulldatasel$season==s)][1])
  USpopls  <- sapply(1:nseas, function(s) USpopdata$USpop[which(USpopdata$agcat==ag & USpopdata$season==s)])
  
  alpha <- c(1,1,1,1) ### Dirichlet prior
  
  data <- list('FSNfluhospnonfatls'=FSNfluhospnonfatls,'FSNfluhospfatls'=FSNfluhospfatls,'FSNpopls'=FSNpopls, 'USpopls' = USpopls,
               'Npideathls'=Npideathls,'oshpideathls'=oshpideathls,
               'nttypearr'=nttypearr,'testposarr'=testposarr,'pcrsens'=pcrsens,
               'lrapidsens'=lrapidsens,'ntotarr'=ntotarr,'alpha'=alpha,'nseas'=nseas)
  
  ptesttot <- t(sapply(1:2, function(k) sapply(1:nseas, function(s) sum(nttypearr[k,,s] +1))))
  
  FSNfluhospls <- FSNfluhospfatls + FSNfluhospnonfatls
  
  pdinit <- FSNfluhospfatls/FSNfluhospls ## initial value of prob. hosp. death
  
  ptestarrinit <- array(0,dim = c(2,4,nseas))
  
  for (t in 1:4) {
    ptestarrinit[,t,] <- (nttypearr[,t,1:nseas] + 1)/ptesttot
  }

  sens1arrinit <- sapply(1:nseas, function(s) c(pcrsens[1],pcrsens[1]))
  logsens2arrinit <- sapply(1:nseas, function(s) c(lrapidsens[1],lrapidsens[1]))
  sens3arrinit <- sapply(1:nseas, function(s) c(.4,.4))
  
  sensarrinit <- array(0,dim = c(2,3,nseas))
  
  for (seas in 1:nseas) {
    for (k in 1:2){
      sensarrinit[k,1,seas] <- sens1arrinit[k,seas]
      sensarrinit[k,2,seas] <- exp(logsens2arrinit[k,seas])
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
      logsens2arr = logsens2arrinit,
      sens3arr = sens3arrinit,
      poshls = poshlsinit,
      USfludeathls = USfludeathlsinit,
      USfluhospls = USfluhosplsinit
    )}
  
  variables <- c('USfludeathls','USfluhospls')
  # variables <- c('sensarr[1:2,1:3,1:7]')
  # variables <- c('ptarr[1:2,1:7]')
  # variables <- c('pt')
  
  setwd(paste0(bfolder,'BEmodels'))
  
  j.model <- jags.model(file=model.file,data=data, inits=inits, n.adapt=nadapt, n.chains=3)
  j.samples<-coda.samples(j.model, variable.names=variables, n.iter=niter, thin = 5) 
  
  cat("\nAge group ",ag," done!!!\n\n")
  
  # summary(j.samples)
  
  codals <- rbind(j.samples[[1]],j.samples[[2]],j.samples[[3]],deparse.level=0)
  
  codaagList[[ag]] <- codals
  codatotls <- codatotls + codals  
}
setwd(paste0(bfolder,'BEmcmc'))

save(codaagList,codatotls,file = outfile)

# sapply(1:nseas,function(s) round(median(codatotls[,s])))
# sapply(1:nseas,function(s) round(median(codatotls[,7 + s])))

#########################################################################################
#########################################################################################
setwd(paste0(bfolder,'BEresults'))

outtablename1 <- "BES hosp beta dirichlet alt 2012-17.txt"
outtablename2 <- "BES mort beta dirichlet alt 2012-17.txt"
#########################################################################################
#########################################################################################
fromyr <- 2012
toyr <- 2017
seasls <- c(sapply(1:nseas, function(seas) paste0(fromyr + (seas-1),'/',substr(fromyr + seas,3,4))),'All Years')

yrls <- fromyr:toyr
firstline <-paste0('Season\t0-4\t5-17\t18-49\t50-64\t65+\tAll Ages')

write.table(firstline,outtablename1,append=T,row.names=FALSE,col.names=FALSE, quote=FALSE)
write.table(firstline,outtablename2,append=T,row.names=FALSE,col.names=FALSE, quote=FALSE)

for (ag in 1:5){
  assign('codaarr',eval(parse(text=paste0('codaagList[[',ag,']]'))))
  
  for (seas in 1:(nseas)){
    assign(paste('eldth',ag,seas,sep = ''), round(as.vector(quantile(codaarr[,seas],prob=c(.5,.05,0.975)))))
  }
  for (seas in 1:(nseas)){
    assign(paste('elhos',ag,seas,sep = ''), round(as.vector(quantile(codaarr[,seas + nseas],prob=c(.5,.05,0.975)))))
  }
}

### summed over all ages
for (seas in 1:nseas){
  assign(paste('eldth',6,seas,sep = ''), round(as.vector(quantile(codatotls[,seas],prob=c(.5,.05,0.975)))))
  assign(paste('elhos',6,seas,sep = ''), round(as.vector(quantile(codatotls[,seas + nseas],prob=c(.5,.05,0.975)))))
}

for (seas in 1:nseas){
  for (ag in 1:6){
    assign(paste('elhos',ag,sep = ''), eval(parse(text = paste0('elhos',ag,seas))))
    assign(paste('eldth',ag,sep = ''), eval(parse(text = paste0('eldth',ag,seas))))
  }
  newline1 <- paste0(seasls[seas],'\t',elhos1[1],' (',elhos1[2],',',elhos1[3],')\t',
                     elhos2[1],' (',elhos2[2],',',elhos2[3],')\t',
                     elhos3[1],' (',elhos3[2],',',elhos3[3],')\t',
                     elhos4[1],' (',elhos4[2],',',elhos4[3],')\t',
                     elhos5[1],' (',elhos5[2],',',elhos5[3],')\t',
                     elhos6[1],' (',elhos6[2],',',elhos6[3],')')
  
  newline2 <- paste0(seasls[seas],'\t',eldth1[1],' (',eldth1[2],',',eldth1[3],')\t',
                     eldth2[1],' (',eldth2[2],',',eldth2[3],')\t',
                     eldth3[1],' (',eldth3[2],',',eldth3[3],')\t',
                     eldth4[1],' (',eldth4[2],',',eldth4[3],')\t',
                     eldth5[1],' (',eldth5[2],',',eldth5[3],')\t',
                     eldth6[1],' (',eldth6[2],',',eldth6[3],')')
  
  setwd(paste0(bfolder,'BEresults'))
  write.table(newline1,outtablename1,append=T,row.names=FALSE,col.names=FALSE, quote=FALSE)
  write.table(newline2,outtablename2,append=T,row.names=FALSE,col.names=FALSE, quote=FALSE)
}

###################################################################################################
###################################################################################################
