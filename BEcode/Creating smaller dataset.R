#########################################################################################
### Aggregating dataset by season, state, ag, testing, test type and test result ########
###  Short version  #####################################################################
#########################################################################################
state <- dataset3$State
seas <- dataset3$season
ag <- dataset3$ag
died <- dataset3$Died
TestedFlu <- dataset3$TestedFlu
TestResult <- dataset3$TestResult

uniqueds2 <- unique(dataset3[,-6])

for (col in 2:length(uniqueds2[1,])){
  uniqueds2[,col] <- as.numeric(as.vector(uniqueds2[,col]))
}
uniqueds2$freq <- uniqueds2$TestResult*0

for (row in seq_along(uniqueds2$State)){
  ds2row <- uniqueds2[row,]
  selind1 <- which(state==ds2row$State & seas==ds2row$season & ag==ds2row$ag & died==ds2row$Died &
                     TestedFlu==ds2row$TestedFlu & TestType==ds2row$TestType & TestResult==ds2row$TestResult)
  uniqueds2$freq[row] <- length(selind1)
  
  agdataset2 <- rbind(agdataset2,row1,deparse.level = 0)
}

for (col in 2:length(agdataset2[1,])){
  agdataset2[,col] <- as.numeric(as.vector(agdataset2[,col]))
}

agdataset2 <- data.frame(agdataset2)

colnames(agdataset2) <- c("state","season","ag","died", "TestedFlu" ,"TestResult","freq")
FluSurvdatashort <- agdataset2
died <- sapply(FluSurvdatashort$died, function(d) ifelse(is.na(d) | d==0,0,1))
FluSurvdatashort$died <- died

FluSurvdatashort <- FluSurvdatashort[order(FluSurvdatashort$state,FluSurvdatashort$season,FluSurvdatashort$ag,
                                 FluSurvdatashort$died,
                                 FluSurvdatashort$TestedFlu,
                                 FluSurvdatashort$TestResult),]
