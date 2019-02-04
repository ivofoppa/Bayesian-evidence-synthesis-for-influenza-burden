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

delind <- which(is.na(dataset$season))
dataset <- dataset[-delind,]

season <- sapply(dataset$season, function(s) switch(which(seasls==s),1,2,3,4,5,6))
dataset$season <- season
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

dataset2 <- dataset[,c("State","season","ag","Died", "TestedFlu" ,"TestType","TestResult","TestType2",
                       "TestResult2","TestType3","TestResult3")]

TestResult <- dataset2$TestResult
TestResult2 <- dataset2$TestResult2

TestResult[which(is.na(TestResult) | TestResult==9 | TestResult==2)] <- 0
TestResult2[which(is.na(TestResult2) | TestResult2==9 | TestResult2==2)] <- 0

TestResult <- sapply(seq_along(TestResult), 
                     function(k) ifelse(TestResult[k]==1,1,
                                        ifelse(TestResult2[k]==1,1,0)))


TestType <- dataset2$TestType
TestType2 <- dataset2$TestType2

TestType[which(is.na(TestType) | TestType==9 | TestType==7)] <- 0
TestType2[which(is.na(TestType2) | TestType2==9 | TestType2==7)] <- 0

TestType <- sapply(seq_along(TestType),function(k) ifelse(TestResult[k]==1,TestType[k],
                                       ifelse(TestResult2[k]==1,TestType2[k],0)))

TestType <- sapply(TestType,function(tt) ifelse(tt>2,3,tt))

dataset2$TestType <- TestType
dataset2$TestResult <- TestResult

dataset3 <- dataset2[,c("State","season","ag","Died", "TestedFlu" ,"TestType","TestResult")]

TestedFlu <- dataset3$TestedFlu

dataset3$TestedFlu[which(is.na(TestedFlu) | TestedFlu>=2)] <- 0
dataset3$TestedFlu[which(TestResult==1)] <- 1

for (col in 2:length(dataset3[1,])){
  dataset3[,col] <- as.numeric(as.vector(dataset3[,col]))
}

#########################################################################################
### Aggregating dataset by season, state, ag, testing, test type and test result ########
#########################################################################################
state <- dataset3$State
seas <- dataset3$season
ag <- dataset3$ag
died <- dataset3$Died
TestedFlu <- dataset3$TestedFlu
TestResult <- dataset3$TestResult
TestType <- dataset3$TestType

uniqueds <- unique(dataset3)

for (col in 2:length(uniqueds[1,])){
  uniqueds[,col] <- as.numeric(as.vector(uniqueds[,col]))
}

agdataset <- NULL

for (k in seq_along(uniqueds$State)){
  dsrow <- uniqueds[k,]
  selind1 <- which(state==dsrow$State & seas==dsrow$season & ag==dsrow$ag & died==dsrow$Died &
                     TestedFlu==dsrow$TestedFlu & TestType==dsrow$TestType & TestResult==dsrow$TestResult)
  row1 <- c(dsrow$State,dsrow$season,dsrow$ag,dsrow$Died,dsrow$TestedFlu,dsrow$TestType,
            dsrow$TestResult,length(selind1))
  
  agdataset <- rbind(agdataset,row1,deparse.level = 0)
}

for (col in 2:length(agdataset[1,])){
  agdataset[,col] <- as.numeric(as.vector(agdataset[,col]))
}

agdataset <- data.frame(agdataset)

colnames(agdataset) <- c("state","season","ag","died", "TestedFlu" ,"TestType","TestResult","freq")
FluSurvdata <- agdataset
died <- sapply(FluSurvdata$died, function(d) ifelse(is.na(d) | d==0,0,1))
FluSurvdata$died <- died

FluSurvdata <- FluSurvdata[order(FluSurvdata$state,FluSurvdata$season,FluSurvdata$ag,
                                 FluSurvdata$died,
                                 FluSurvdata$TestedFlu,
                                 FluSurvdata$TestType,
                                 FluSurvdata$TestResult),]

# FluSurvdatashort <- FluSurvdata[,c('state','season','ag','died','TestedFlu','TestResult','freq')]
# FluSurvdatashort <- FluSurvdatashort[order(FluSurvdatashort$state,FluSurvdatashort$season,FluSurvdatashort$ag,
#                                  FluSurvdatashort$died,
#                                  FluSurvdatashort$TestedFlu,
#                                  FluSurvdatashort$TestResult),]

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
