#########################################################################################
###  Code to sum NCHS data over state, season, age group, cod, place of death  ##########
#########################################################################################
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,86),'65+'))

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
  s <- which(seasls==seas)
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

seas <- sapply(oshdat$Season, function(s) which(seasls==s))
oshdat$seas <- seas

ag <- sapply(oshdat$Age, function(ag) which(agels==ag))
oshdat$ag <- ag

oshdat$freq <- oshdat$COUNT
oshdat$state <- oshdat$State

oshdat <- oshdat[,c('seas','state','ag','osh','freq')]
#########################################################################################
### Age specific estimates from Millman et al. EID 2015, 21 (9):
cipcrlist <- list(c(95.0,82,98.7),c(95.0,82,98.7),c(94.1,81.1,98.7),c(94.1,81.1,98.7),c(86.1,79.6,92.7))
cirapidlist <- list(c(66.7,61.3,71.7),c(66.7,61.3,71.7),c(53.9,47.8,59.8),c(53.9,47.8,59.8),c(20.1,8.8,41.4))
sensdata <- list(cipcrlist,cirapidlist)
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
## "0-4y"   "5-17y"  "18-49y" "50-64y" "65+y"
ag <- sapply(dataset$Age, function(ag) ifelse(ag %in% c(0:4),1,
                                              ifelse(ag %in% c(5:17),2,
                                                     ifelse(ag %in% c(18:49),3,
                                                            ifelse(ag %in% c(50:64),4,
                                                                   ifelse(ag %in% c(65:110),5,NA))))))
dataset$ag <- ag
dataset <- dataset[which(!is.na(ag)),]
dataset$Died <- sapply(dataset$Died,function(d) ifelse(is.na(d),0,d))
#########################################################################################
### Organizing dataset by season, state, ag, testing, test type and test result #########
#########################################################################################
# "State", "TestedFlu" ,"TestType","TestResult","TestType2","TestResult2","TestType3","TestResult3",
#"Died","season","othtest","ag"
dataset2 <- dataset[,c("State","season","ag", "TestedFlu" ,"TestType","TestResult","TestType2","TestResult2","TestType3","TestResult3",
                       "Died")]
TestType <- sapply(seq_along(dataset2$TestType), function(t) ifelse(dataset2$TestResult[t]==1,
                                                                             dataset2$TestType[t],
                                                                             ifelse(dataset2$TestResult2[t]==1,
                                                                                    dataset2$TestType2[t],
                                                                                    ifelse(dataset2$TestResult3[t]==1,
                                                                                           dataset2$TestType3[t],
                                                                                           dataset2$TestType[t]))))

TestResult <- sapply(seq_along(dataset2$TestResult), function(t) ifelse(dataset2$TestResult[t]==1,
                                                                             dataset2$TestResult[t],
                                                                             ifelse(dataset2$TestResult2[t]==1,
                                                                                    dataset2$TestResult2[t],
                                                                                    ifelse(dataset2$TestResult3[t]==1,
                                                                                           dataset2$TestResult3[t],
                                                                                           dataset2$TestResult[t]))))
TestType <- sapply(TestType,function(tt) ifelse(tt==1,1,ifelse(tt==2,2,ifelse(tt>2,3,0))))
dataset2$TestType <- TestType
dataset2$TestResult <- TestResult

dataset2 <- dataset[,c("State","season","ag", "TestedFlu" ,"TestType","TestResult","Died")]

dataset2$TestedFlu[which(is.na(dataset2$TestedFlu) | dataset2$TestedFlu==2)] <- 0
dataset2$TestedFlu[which(dataset2$TestResult==1)] <- 1

dataset2$TestResult[which(is.na(dataset2$TestResult) | dataset2$TestResult==2)] <- 0
#########################################################################################
### Aggregating dataset by season, state, ag, testing, test type and test result ########
#########################################################################################
State <- dataset2$State
seas <- dataset2$seas
ag <- dataset2$ag
TestedFlu <- dataset2$TestedFlu
TestedType <- dataset2$TestedType
TestResult <- dataset2$TestResult
Died <- dataset2$Died

agdataset <- NULL

for (st in statels){
  for (s in 1:6){
    for (a in 1:5){
      for (died in 0:1){
        for (t in 0:1){
          if (t==1){
            for (tt in 1:3){
              for (res in 0:1){
                selind <- which(State==st & seas==s & ag==a & Died==died & TestedFlu==1 & TestType==tt & TestResult==res)
                row <- c(st,s,a,died,t,tt,res,length(selind))
              }
            }
          } else {
            selind <- which(State==st & seas==s & ag==a & Died==died & TestedFlu==0)
            row <- c(st,s,a,died,0,NA,NA,length(selind))
          }
          agdataset <- rbind(agdataset,row,deparse.level = 0)
        }
      }
      
    }
  }
}

for (col in 2:length(agdataset[1,])){
  agdataset[,col] <- as.numeric(agdataset[,col])
}

agdataset <- data.frame(agdataset)

colnames(agdataset) <- c("State","season","ag", "TestedFlu" ,"TestType","TestResult","Died","freq")
FluSurvdata <- agdataset
#########################################################################################
###  Reading-in denominator data ########################################################
#########################################################################################
library(readxl)
#,
statels2 <- c("CA","GA","OR","CT","CO","NM","NYA","NYR","MN","TN","MI")

setwd(paste0(bfolder,'BEdata'))
agseaspop <- NULL

for (s in 1:7){
  seas <- seaslabls[[s]]
  xlsFile <- paste0('NCHS ',seas,' population estimates.xls')
  popfile <- read_excel(xlsFile, sheet = 1)
  for (st in statels){
    stpopselind <- which(substr(statels2,1,2)==st)+1
    for (ag in 1:5){
      aglim <- agecatlist[[ag]][[1]]
      drow <- c(s,ag,st,sum(popfile[aglim,stpopselind]))
      agseaspop <- rbind(agseaspop,drow,deparse.level = 0)
    }
  }
}
colnames(agseaspop) <- c('seas','ag','state','pop')
agseaspop <- data.frame(agseaspop)
agseaspop <- agseaspop[,c('state','seas','ag','pop')]
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
outfname <- 'FluSURV-NET_burden.RData'
fulldata <- list(FluSurvdata=FluSurvdata,sensdata=sensdata,agseaspop=agseaspop,oshdat=oshdat)

save(fulldata,file = outfname)
#########################################################################################
#########################################################################################
#########################################################################################
