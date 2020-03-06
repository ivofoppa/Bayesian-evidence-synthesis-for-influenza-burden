library(readxl)
library(lubridate) ## for extracting months etc. from dates
#########################################################################################
#########################################################################################
rm(list = ls())
bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

nseas <- 8
seasvec <- sapply(1:nseas,function(y) paste0(y + 9,y + 10)) ### for reading-in NCHS data
#########################################################################################
###  Processing FIPS county data ########################################################
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
FIPSdata <- read.csv("county_FIPS.csv")

#########################################################################################
### Collecting FSN denominator data #####################################################
#########################################################################################
yearCountyArr <- array(,dim = c(0,3))

for (seas in 1:nseas) {
  # fipsls <- NULL
  fname <- paste0("NCHS ",seasvec[seas] ," population estimates.xls")
  dataset <- read_excel(fname,sheet = 'County')
  statevec <- unique(substr(dataset$state,1,2))
  for (st in statevec) {
    stdataset <- dataset[which(substr(dataset$state,1,2)==st),]
    countyvec <- stdataset$county
    for (cnty in countyvec) {
      selind <- which(FIPSdata$Name==cnty & FIPSdata$State==st)
      fipscd <- FIPSdata$FIPS[selind]
      fipscd <- ifelse(fipscd > 9999,toString(fipscd),paste0(0,toString(fipscd)))
      fips <- as.numeric(substr(fipscd,3,5))
      yearCountyArr <- rbind(yearCountyArr,c(seas,st,fips),deparse.level = 0)
      # fipsls <- c(fipsls,fips)
    }
  }
}
colnames(yearCountyArr) <- c('Season','State','FIPS')

FIPScodes <- data.frame(yearCountyArr)

#########################################################################################
#########################################################################################
setwd(paste0(bfolder,'BEdata'))
outfname <- "FIPScodeSeas.txt"
write.table(FIPScodes,file = outfname,row.names = F,quote = F)
