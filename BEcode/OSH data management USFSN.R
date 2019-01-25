#########################################################################################
###  Code to sum NCHS data over state, season, age group, cod, place of death  ##########
#########################################################################################
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

#########################################################################################
#########################################################################################
### Manipulating osh data from FluSurv-NET
#########################################################################################
#########################################################################################
infname <- 'Deaths_Ivo_20190124.csv'
setwd(paste0(bfolder,'BEdata'))

oshdat <- read.csv(infname)

for (col in seq_along(oshdat[1,])){
  oshdat[,col] <- unlist(as.vector(oshdat[,col]))
}

oshdat$osh <- sapply(oshdat$Death, function(dth) ifelse(dth=='Died 30 days post-discharge',1,0))
seasls <- unique(oshdat$Season)
statels <- unique(oshdat$State)
agels <- unique(oshdat$Age)
oshls <- unique(oshdat$Season)
#########################################################################################
##  to combine data from the two NY sites ###############################################
#########################################################################################
for (seas in seasls){
  for (ag in agels){
    for (osh in 0:1){
      selind1 <- which(oshdat$Season==seas & oshdat$Age==ag & oshdat$osh==osh & oshdat$State=='NYA')
      selind2 <- which(oshdat$Season==seas & oshdat$Age==ag & oshdat$osh==osh & oshdat$State=='NYR')
      row1 <- oshdat[selind1,]
      row2 <- oshdat[selind2,]
      rowcomb <- row1
      rowcomb[5] <- row1[5] + row2[5]
      rowcomb[2] <- "NY"
      oshdat <- oshdat[-c(selind1,selind2),]
      oshdat <- rbind(oshdat,rowcomb)
      }
  }
}
#########################################################################################

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
### For creation of a season variable in the FluSurv-NET data set #######################
mmwrstrt <- as.vector(mmwrdat$mmwrstrt)
mmwrstrt <- as.Date(mmwrstrt,"%m/%d/%Y")

mmwrend <- as.vector(mmwrdat$mmwrend)
mmwrend <- as.Date(mmwrend,"%m/%d/%Y")

for (seas in 1:6){
  ls <- which(DateHosp >= mmwrstrt[seaslist[[seas]][1]] & DateHosp <= mmwrend[seaslist[[seas]][2]])
  dataset$season[ls] <- seas
}

#########################################################################################

statels <- unique(dataset$State)
seaslabls <- sapply(1:7, function(s) paste0(10 + s-1,10 + s))
#########################################################################################
###  Reading-in denominator data ########################################################
#########################################################################################
library(readxl)
# 
for (s in 1:7){
  seas <- seaslabls[[s]]
  xlsFile <- paste0('NCHS ',seas,' population estimates.xls')
  popfile <- read_excel(xlsFile, sheet = 1)
  ls <- names(popfile)
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
