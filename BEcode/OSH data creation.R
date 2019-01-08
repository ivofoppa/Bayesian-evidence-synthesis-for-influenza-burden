#########################################################################################
###  Code to sum NCHS data over state, season, age group, cod, place of death  ##########
#########################################################################################
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

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

agcat <- 4

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
setwd(paste0(bfolder,'BEdata'))
mmwrdat0<-read.table("mmwrweeks.txt",header=T)
mmwrdat<-mmwrdat0[which((mmwrdat0$wyear==fromyr&mmwrdat0$week>=fromwk)|(mmwrdat0$wyear>fromyr&mmwrdat0$wyear<toyr)|(mmwrdat0$wyear==toyr&mmwrdat0$week<=towk)),]

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
dataset$DateHosp <- DateHosp

statels <- unique(dataset$State)
seaslabls <- sapply(1:7, function(s) paste0(10 + s-1,10 + s))
#########################################################################################
### NCHS data set: Weekly mortality by age group, week, diagnostic group and ############
### place of death for proportion deaths outside the hospital                ############
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
dataset2 <- data.frame(read.csv('mort2010_16_week_5agstate.csv'))
dataset2$oshospdth <- sapply(dataset2$placdth, function(plce) ifelse(plce==1,0,1) )

dataset2$DVD <- sapply(dataset2$studywk, function(week) mmwrdat0$dvdweek[which(mmwrdat0$studywk==week)])
### Changing vars into simple values
for (k in 1:length(dataset2[1,])) {
  dataset2[,k] <- as.vector(dataset2[,k])
}
stateselind <- which(dataset2$stateoc %in% statels & dataset2$DVD %in% DVDun) ## restrict data to states represented in FluSurv-NET data
dataset2 <- dataset2[stateselind,]
###  Creating a season variable #########################################################
seasvar <- NULL
for (seas in 1:nseas){
  seaslim <- seaslist[[seas]]
  seasselind <- which(dataset2$DVD %in% DVDun[seaslim[1]:seaslim[2]])
  seasvar[seasselind] <- seas
}

dataset2$season <- seasvar

###  summing over seasons ###############################################################
dataset2cum <- NULL
for (seas in 1:nseas){
  for (ag in 1:5){
    for (state in statels){
      for (osh in 0:1){
        selind <- which (dataset2$season==seas & dataset2$agecat==ag & dataset2$stateoc==state & dataset2$oshospdth==osh)
        dataset2sel <- dataset2[selind,]
        cumrow <- c(seas,state,ag,sum(dataset2sel$rcu),sum(dataset2sel$pi),osh)
        dataset2cum <- rbind(dataset2cum, cumrow, deparse.level = 0)
      }
    }
  }
}

for (col in c(1,3:6)){
  dataset2cum[,col] <- as.numeric(dataset2cum[,col])
}
colnames(dataset2cum) <- c('season','state','ag','rcu','pi','osh')
dataset2cum <- data.frame(dataset2cum)
#########################################################################################
###  Reading-in denominator data ########################################################
#########################################################################################
library(readxl)
# 
dataset2cum$pop <- dataset2cum$rcu*0


for (s in 1:7){
  seas <- seaslabls[[s]]
  xlsFile <- paste0('NCHS ',seas,' population estimates.xls')
  popfile <- read_excel(xlsFile, sheet = 1)
  for (st in statels){
    for (ag in 1:5){
      aglim <- agecatlist[[ag]][[1]]
      stpopselind <- which(substr(ls,1,2)==st)
      if (length(stpopselind) > 1){
        stpopdat <- cbind(popfile$age, rowSums(popfile[,stpopselind]),deparse.level = 0)
      } else stpopdat <- cbind(popfile$age, popfile[,stpopselind],deparse.level = 0)
      agpopselind <- which(stpopdat[,1] >= aglim[1] & stpopdat[,1] <= aglim[2])
      agpop <- sum(stpopdat[agpopselind,2])
      
    }
  }
  
}

#########################################################################################
setwd(paste0(bfolder,'BEdata'))
outfname <- 'OSHdeath_2010-2016,FluSURV-NET.csv'
write.csv(dataset2cum,file = outfname,quote = F)
#########################################################################################
#########################################################################################
#########################################################################################
