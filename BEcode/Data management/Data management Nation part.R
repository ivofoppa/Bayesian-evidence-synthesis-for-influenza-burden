#########################################################################################
###  Code to sum NCHS data over state, season, age group, cod, place of death  ##########
#########################################################################################
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,110),'65+'))

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

# agcat <- 4
# 
# agecatls <- agecatlist[[agcat]]
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

DVD <- mmwrdat$dvdweek

N <- length(mmwrdat$dvdweek)
seasbeg <- which(mmwrdat$week==fromwk)
seasend <- which(mmwrdat$week==towk)
time <- (1:N)/N

seaslist <- list()
for (seas in 1:nseas){
  seaslist[[seas]] <- c(DVD[seasbeg[seas]],DVD[seasend[seas]])
}
#########################################################################################
###  FluSurv-NET testing data set #######################################################
#########################################################################################
### Data managment for detection multiplier
setwd(paste0(bfolder,'BEdata'))

dataset <- data.frame(read.csv('20181214_1415burden_bysite_simple.csv'))

selind <- which(dataset$State%in%c("CA","CO","NM","NY","OR"))
datasetpart <- dataset[selind,]

### Changing vars into simple values
for (k in 4:length(datasetpart[1,])) {
  datasetpart[,k] <- as.integer(as.vector(datasetpart[,k]))
}

DateHosp <- as.Date(datasetpart$DateHosp,"%m/%d/%Y")
datasetpart$DateHosp <- DateHosp

mmwrdat0$mmwrstrt <- as.Date(mmwrdat0$mmwrstrt,"%m/%d/%Y")
mmwrdat0$mmwrend <- as.Date(mmwrdat0$mmwrend,"%m/%d/%Y")

DVD <- sapply(DateHosp, function(d) mmwrdat0$dvdweek[which(mmwrdat0$mmwrstrt<=d & mmwrdat0$mmwrend >=d)])
datasetpart$DVD <- DVD
### selecting relevant time period
selind <- which(datasetpart$DVD>=as.numeric(paste0(fromyr,fromwk)) & datasetpart$DVD <= as.numeric(paste0(toyr,towk)))
datasetpart <- datasetpart[selind,]

seaslabls <- sapply(1:6, function(s) paste0(10 + s-1,10 + s))

### modified season: 1,...,6
season2 <- datasetpart$season
for (k in 1:6) {
  selind <- which(datasetpart$DVD >= seaslist[[k]][1] & datasetpart$DVD <= seaslist[[k]][2])
  season2[selind] <- k
}
datasetpart$season <- season2
### Defining agecat
agecat <- datasetpart$Age
for (ag in 1:5) {
  selind <- which(datasetpart$Age %in% c(agecatlist[[ag]][[1]][1]:agecatlist[[ag]][[1]][2]))
  agecat[selind] <- round(ag)
}
datasetpart$agecat <- agecat

TestTyperec <- sapply(datasetpart$TestType, function(t) ifelse(t==1,1,ifelse(t==2,2,ifelse(t%in%3:6,3,0))))
datasetpart$TestTyperec <- as.vector(TestTyperec)

datasetcum <- NULL

for (seas in 1:6) {
  for (ag in 1:5) {
 ## died
    selseasagind1 <- which(datasetpart$agecat==ag & 
                             datasetpart$season==seas & datasetpart$Died==1)
    datasetagseas1 <- datasetpart[selseasagind1,]
    if (length(selseasagind1) > 0) {
      for (tt in 1:3) {
        numtest1 <- length(which(datasetagseas1$TestedFlu==1 &
                                   datasetagseas1$TestTyperec==tt & datasetagseas1$TestResult==1))
        numtest0 <- length(which(datasetagseas1$TestedFlu==1 &
                                   datasetagseas1$TestTyperec==tt & datasetagseas1$TestResult!=1))
        
        el1 <- c(seas,ag,1,1,tt,1,numtest1)
        el0 <- c(seas,ag,1,1,tt,0,numtest0)
        datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
      }
    } else {
      el1 <- c(seas,ag,1,1,tt,1,0)
      el0 <- c(seas,ag,1,1,tt,0,0)
      datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
    }
    num <- length(which(datasetagseas1$TestedFlu!=1))
    el <- c(seas,ag,1,0,0,0,num)
    datasetcum <- rbind(datasetcum,el,deparse.level = 0)
    ##did not die
    selseasagind0 <- which(datasetpart$agecat==ag & 
                             datasetpart$season==seas & datasetpart$Died!=1)
    datasetagseas0 <- datasetpart[selseasagind0,]
    if (length(selseasagind1) > 0) {
      for (tt in 1:3) {
        numtest1 <- length(which(datasetagseas0$TestedFlu==1 &
                                   datasetagseas0$TestTyperec==tt & datasetagseas0$TestResult==1))
        numtest0 <- length(which(datasetagseas0$TestedFlu==1 & 
                                   datasetagseas0$TestTyperec==tt & datasetagseas0$TestResult!=1))
        
        el1 <- c(seas,ag,0,1,tt,1,numtest1)
        el0 <- c(seas,ag,0,1,tt,0,numtest0)
        datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
      } 
    } else {
        el1 <- c(seas,ag,0,1,tt,1,0)
        el0 <- c(seas,ag,0,1,tt,0,0)
        datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
    }
    
    num <- length(which(datasetagseas0$TestedFlu!=1))
    el <- c(seas,ag,0,0,0,0,num)
    datasetcum <- rbind(datasetcum,el,deparse.level = 0)
  }
}
colnames(datasetcum) <- c('season','agecat','died','TestedFlu','TestType','TestResult','freq')
FSNcumtestdata <- data.frame(datasetcum)
#########################################################################################
### FSN data by outcome (fatal, non-fatal)   ############################################
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
FSNfname <- "foppa_14MAY19.csv"
FSNdata <- read.csv(FSNfname)
agecatls <- unique(FSNdata$Age)
agls <- as.vector(agecatls)

for (ag in agecatls) {
  selind <- which(FSNdata$Age==ag)
  agls[selind] <- which(agecatls==ag)
}
FSNdata$agecat <- as.numeric(agls)

occatls <- unique(as.vector(FSNdata$Outcome))

ocls <- as.vector(FSNdata$Outcome)

for (oc in occatls) {
  selind <- which(FSNdata$Outcome==oc)
  ocls[selind] <- which(occatls==oc)
}

FSNdata$outcome <- 2 - as.numeric(ocls)

seasls <- unique(FSNdata$Season)

for (seas in seasls) {
  selind <- which(FSNdata$Season==seas)
  FSNdata$Season[selind] <- which(seasls==seas)
}
FSNdata <- FSNdata[,c("Season","State","agecat","outcome","Count")]
FSNdata$state <- as.vector(FSNdata$state)
colnames(FSNdata) <- c("season","state","agecat","died","freq")

agls <- unique(FSNdata$agecat)
dls <- unique(FSNdata$died)
seasls <- unique(FSNdata$season)
statels <- unique(FSNdata$state)

dataarr <- NULL
for (seas in seasls) {
  for (ag in agls) {
    for (d in dls) {
      selind <- which(FSNdata$season==seas & FSNdata$agecat==ag & FSNdata$died==d)
      row <- c(season=seas,agecat=ag,died=d,freq=sum(FSNdata$freq[selind]))
      dataarr <- rbind(dataarr,row,deparse.level = 0)
    }
  }
}

FSNcumdata <- data.frame(dataarr)
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
###  Creating a season variable #########################################################
DVDun <- unique(DVD)
seasvar <- dataset2$DVD*0
for (seas in 1:nseas){
  seaslim <- seaslist[[seas]]
  seasselind <- which(dataset2$DVD %in% c(seaslim[1]:seaslim[2]))
  seasvar[seasselind] <- seas
}

dataset2$season <- seasvar

###  summing over seasons ###############################################################
dataset2cum <- NULL
for (seas in 1:nseas){
  for (ag in 1:5){
    for (osh in 0:1){
      selind <- which (dataset2$season==seas & dataset2$agecat==ag & dataset2$oshospdth==osh)
      dataset2sel <- dataset2[selind,]
      cumrow <- c(seas,ag,sum(dataset2sel$rcu),sum(dataset2sel$pi),osh)
      dataset2cum <- rbind(dataset2cum, cumrow, deparse.level = 0)
    }
  }
}

for (col in c(1:5)){
  dataset2cum[,col] <- as.numeric(dataset2cum[,col])
}
colnames(dataset2cum) <- c('season','agecat','rcu','pi','osh')
OSHcumdata <- data.frame(dataset2cum)
#########################################################################################
###  Reading-in denominator data ########################################################
#########################################################################################
# 

popdata <- read.csv("FSNpop.csv")

#########################################################################################
setwd(paste0(bfolder,'BEdata'))
outfname <- 'FluSURV-NET-Nation_part.Rdata'
save(FSNcumtestdata,FSNcumdata,OSHcumdata,popdata,sensdata,file = outfname)
#########################################################################################
#########################################################################################
# #########################################################################################
# bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'
# setwd(paste0(bfolder,'BEdata'))
# infname <- 'FluSURV-NET.Rdata'
# load(infname)
