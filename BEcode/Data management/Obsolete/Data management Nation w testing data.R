library(readxl)
library(lubridate) ## for extracting months etc. from dates
#########################################################################################
#########################################################################################
rm(list = ls())
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

nseas <- 8

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

agcat <- 5
agernge <- agecatlist[[agcat]][[1]]
agernge[2] <- min(agernge[2],85) 
agecatlist[[agcat]][[1]] <- agernge

agelims <- sapply(1:5,function(ag) agecatlist[[ag]][[1]])

seasvec <- sapply(1:nseas,function(y) paste0(y + 9,y + 10)) ### for reading-in NCHS data
### Load data 
setwd(paste0(bfolder,'BEdata'))
#########################################################################################
###  Processing testing data ############################################################
#########################################################################################
Testdata <- read.csv("ivo_burden_20190523.csv")
myvars <- c("State","DateHosp","TestedFlu","TestResult2","DateDeath","Gender","Test","pcr","DateHosp","TestType",
             "TestType3","HospID","ICU","season","rapid","othtest","Age","TestResult","TestResult3","Facility","dob",
             "TestType2","Died")

Testdata <- Testdata[myvars]
for (col in 1:length(myvars)) {
  Testdata[,col] <- as.vector(Testdata[,col])
}

hospdates <- Testdata$DateHosp
hospdates <- as.Date(hospdates,"%d%B%Y")

### Deleting obs with missing hosp date
delind <- which(hospdates < "2010-10-01" | hospdates > "2018-06-01" | is.na(hospdates))
Testdata <- Testdata[-delind,]

hospdates <- Testdata$DateHosp
hospdates <- as.Date(hospdates,"%d%B%Y")

months <- month(as.POSIXlt(hospdates, format="%Y-%m-%Y"))
years <- year(as.POSIXlt(hospdates, format="%Y-%m-%Y"))

for (k in seq_along(hospdates)) {
  mo <-  months[k]
  yr <- as.numeric(substr(years[k],3,4))
  seas <- ifelse(mo > 6,paste0(yr,yr+1),paste0(yr-1,yr))
  Testdata$season[k] <- which(seasvec==seas)
}

age <- Testdata$Age
agcatvec <- as.numeric(sapply(age,function(a) which(agelims[1,] <= a & agelims[2,] >= a )))
Testdata$agcat <- agcatvec
Testdata$Died <- as.numeric(sapply(Testdata$Died, function(d) ifelse(d==1,1,0)))
## Aggregating
## variables to create:
## separately by died
## # tested by test type (PCR,AG,other,none)
## # results per test type
TestdataCum <- NULL
for (seas in 1:nseas) {
  for (ag in 1:5) {
    for (d in 1:0) {
      selind <- which(Testdata$season==seas & Testdata$agcat==ag & Testdata$Died==d)
      dat <- Testdata[selind,]
      ntest1 <- length(which(dat$TestType==1))
      ntest2 <- length(which(dat$TestType==2))
      ntest3 <- length(which(dat$TestType>2))
      ntest0 <- length(which(is.na(dat$TestType)))
      
      testres1 <- length(which(dat$TestType==1 & dat$TestResult==1))
      testres2 <- length(which(dat$TestType==2 & dat$TestResult==1))
      testres3 <- length(which(dat$TestType>2 & dat$TestResult==1))
      
      row <- c(seas,ag,d,ntest1,ntest2,ntest3,ntest0,testres1,testres2,testres3)
      TestdataCum <- rbind(TestdataCum,row,deparse.level = 0)      
    }
  }
}
TestdataCum <- data.frame(TestdataCum)
colnames(TestdataCum) <- c("season","agcat","died","ntest1","ntest2","ntest3","ntest0","testres1","testres2","testres3")

#########################################################################################
### FSN data by outcome (fatal, non-fatal)   ############################################
#########################################################################################
FSNfname <- "foppa_14MAY19.csv"
FSNdata <- read.csv(FSNfname)

agecatls <- unique(FSNdata$Age)
agls <- as.vector(agecatls)

for (ag in agecatls) {
  selind <- which(FSNdata$Age==ag)
  agls[selind] <- which(agecatls==ag)
}
FSNdata$agcat <- as.numeric(as.vector(agls))

occatls <- unique(as.vector(FSNdata$Outcome))
ocls <- as.vector(FSNdata$Outcome)

for (oc in occatls) {
  selind <- which(FSNdata$Outcome==oc)
  ocls[selind] <- which(occatls==oc)
}

FSNdata$outcome <- 2 - as.numeric(as.vector(ocls))

seasls <- unique(FSNdata$Season)

for (seas in seasls) {
  selind <- which(FSNdata$Season==seas)
  FSNdata$Season[selind] <- which(seasls==seas)
}

FSNdata <- FSNdata[,c("Season","State","agcat","outcome","Count")]

colnames(FSNdata) <- c("season","state","agcat","died","freq")
### Aggregating over states
FSNdataCum <- NULL
for (seas in 1:nseas) {
  for (ag in 1:5) {
    for (d in 0:1) {
      selind <- which(FSNdata$season==seas & FSNdata$agcat==ag & FSNdata$died==d)
      num <- sum(FSNdata$freq[selind])
      FSNdataCum <- rbind(FSNdataCum,c(seas,ag,d,num),deparse.level = 0)
    }
  }
}
FSNdataCum <- data.frame(FSNdataCum)
colnames(FSNdataCum) <- c("season","agcat","died","freq")
#########################################################################################
### US population data (NCHS)  ##########################################################
#########################################################################################
popfname <- "census_NCHS.xls"
poptotdata <- NULL

yearls <- 2010:2017

for (k in seq_along(yearls)) {
  yrstr <- paste0("_",yearls[k])
  popdata0 <- read_excel(popfname,sheet = yrstr)
  popvec <- as.vector(unlist(popdata0[,2]))
  agvec <- 1:5
  seasls <- rep(k,5)
  subdata <- cbind(seasls,agvec,popvec,deparse.level = 0)
  poptotdata <- rbind(poptotdata,subdata,deparse.level = 0)
}
USpoptotdata <- data.frame(poptotdata)
colnames(USpoptotdata) <- c("season","agcat","USpop")
#########################################################################################
### Collecting FSN denominator data #####################################################
#########################################################################################
FSNpopdata <- NULL
for (seas in 1:nseas) {
  fname <- paste0("NCHS ",seasvec[seas] ," population estimates.xls")
  dataset <- read_excel(fname)
  for (ag in 1:5) {
    agel <- agelims[,ag]
    FSNpop <- sum(dataset[(agel[1]:agel[2]) + 1,]) ## sums over states/age ranges of relevance 
    row <- c(seas,ag,FSNpop)
    FSNpopdata <- rbind(FSNpopdata,row,deparse.level = 0)
  }
}
FSNpopdata <- data.frame(FSNpopdata)
colnames(FSNpopdata) <- c("season","agcat","FSNpop")
### example for season 4
#########################################################################################
### Creating final data set  ############################################################
#########################################################################################
#########################################################################################
Fulldata <- NULL
for (seas in 1:nseas) {
  for (ag in 1:5) {
    Fpopselind <- which(FSNpopdata$seas==seas & FSNpopdata$ag==ag )
    USpopselind <- which(USpoptotdata$seas==seas & USpoptotdata$ag==ag )
    Fpop <- FSNpopdata[Fpopselind,"FSNpop"]
    USpop <- USpoptotdata[USpopselind,"USpop"]
    
    for (d in 0:1) {
      Testselind <- which(TestdataCum$season==seas & TestdataCum$agcat==ag&TestdataCum$died==d)
      FSNselind <- which(FSNdataCum$season==seas & FSNdataCum$agcat==ag&FSNdataCum$died==d)
      
      tdat <- TestdataCum[Testselind,c("ntest1","ntest2","ntest3","ntest0","testres1","testres2","testres3")]
      FSNoutcome <- FSNdataCum[FSNselind,"freq"]
      
      row <- unlist(c(seas,ag,d,tdat,FSNoutcome,Fpop,USpop))
      Fulldata <- rbind(Fulldata,row,deparse.level = 0)
    }
  }
}

Fulldata <- data.frame(Fulldata)
colnames(Fulldata) <- c("season","agcat","died","ntest1","ntest2","ntest3","ntest0","testres1",
                        "testres2","testres3","Noutcome","FSNpop","USpop")  

#########################################################################################
### Reading-in OSH data      ############################################################
#########################################################################################
OSHfname <- "mort2010_17_season_osh.csv"
OSHdata <- read.csv(OSHfname)
#########################################################################################
nseas <- 7

cipcrlist <- list(c(95.0,82,98.7),c(95.0,82,98.7),c(94.1,81.1,98.7),c(94.1,81.1,98.7),c(86.1,79.6,92.7))
cirapidlist <- list(c(66.7,61.3,71.7),c(66.7,61.3,71.7),c(53.9,47.8,59.8),c(53.9,47.8,59.8),c(20.1,8.8,41.4))
sensdata <- list(cipcrlist,cirapidlist)
#########################################################################################
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
### Load data 
infname <- 'FluSURV-NET-Nation-osh-7seas.RData'

save(Fulldata,OSHdata,nseas,cipcrlist,sensdata,file = infname)
