library(readxl)
#########################################################################################
#########################################################################################
rm(list = ls())
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'
### Load data 
setwd(paste0(bfolder,'BEdata'))

fname <- paste0('AgeseasoncodaList_Hosp.RData')
load(fname)

USpopdata <- read.csv("USpop.csv")
USpopdata$age <- as.vector(USpopdata$age)

delind <- which(as.factor(USpopdata$age)=="")
USpopdata <- USpopdata[-delind,]

selind <- which(as.factor(USpopdata$age)=="85+")
USpopdata$age[selind] <- "85"

USpopdata$age <- as.numeric(USpopdata$age)

agecatlist <- list(list(c(0,4),'<5'),
                   list(c(5,17),'5 to 17'),
                   list(c(18,49),'18-49'),
                   list(c(50,64),'50-64'),
                   list(c(65,120),'65+'))

seasvec <- sapply(1:nseas,function(y) paste0(y + 9,y + 10)) ### for reading-in NCHS data

estarr <- array(0,dim = c(5,nseas))
ageestList <- list()

for (agcat in 1:5) {
  agernge <- as.numeric(agecatlist[[agcat]][[1]])
  agernge[2] <- min(agernge[2],85) 
  
 estarr <- NULL
  seasoncodaList <- AgeseasoncodaList[[agcat]]
  
  for (k in 1:nseas) {
    fname <- paste0("NCHS ",seasvec[k] ," population estimates.xls")
    dataset <- read_excel(fname)
    
    stateredvec <- seasoncodaList[[k]][["states"]]
    colnmshort <- as.vector(sapply(colnames(dataset), function(nm) substr(nm,1,2)))
    colselind <- which(colnmshort%in%stateredvec)
    FSNpop <- sum(dataset[(agernge[1]:agernge[2]) + 1,colselind]) ## sums over states/age ranges of relevance 
    # FSNpopvec <- c(FSNpopvec,FSNpop)  
    # 
    year <- eval(parse(text = paste0("20",substr(seasvec[k],1,2))))
    USpop <- sum(USpopdata$pop[which(USpopdata$year==year & (USpopdata$age>=agernge[1] & USpopdata$age<=agernge[2]))])
    # USpopvec <- c(USpopvec,USpop)
    
    FSMest <- quantile(seasoncodaList[[k]][[2]],probs = c(.025,.5,.975))
    
    estarr <- rbind(estarr,round(FSMest/FSNpop*USpop))
  }
  ageestList[[agcat]] <- estarr
}
#########################################################################################
#########################################################################################

pd <- c(0.155830534,0.142431389,0.207497822,0.24405856,0.335799256,0.31614653)

