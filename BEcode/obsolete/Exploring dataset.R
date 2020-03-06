#########################################################################################
### Aggregating dataset by season, state, ag, testing, test type and test result ########
###  Short version, no test type, for checking  #########################################
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

State <- as.vector(uniqueds2$State)
uniqueds2$State <- State

uniqueds2$freq <- uniqueds2$TestResult*0

for (row in seq_along(uniqueds2$State)){
  ds2row <- uniqueds2[row,]
  selind1 <- which(state==ds2row$State & seas==ds2row$season & ag==ds2row$ag & died==ds2row$Died &
                     TestedFlu==ds2row$TestedFlu & TestResult==ds2row$TestResult)
  uniqueds2$freq[row] <- length(selind1)
}

for (col in 2:length(uniqueds2[1,])){
  uniqueds2[,col] <- as.numeric(as.vector(uniqueds2[,col]))
}

colnames(uniqueds2) <- c("state","season","ag","died", "TestedFlu" ,"TestResult","freq")
FluSurvdatashort <- uniqueds2
died <- sapply(FluSurvdatashort$died, function(d) ifelse(is.na(d) | d==0,0,1))
FluSurvdatashort$died <- died

FluSurvdatashort <- FluSurvdatashort[order(FluSurvdatashort$state,FluSurvdatashort$season,FluSurvdatashort$ag,
                                 FluSurvdatashort$died,
                                 FluSurvdatashort$TestedFlu,
                                 FluSurvdatashort$TestResult),]


uniqueds2b <- uniqueds2[which(uniqueds2$TestedFlu==1),-5] 

uniqueds3 <- unique(uniqueds2b[,c(1:4)])
uniqueds3 <- uniqueds3[order(uniqueds3$state,uniqueds3$season,uniqueds3$ag,uniqueds3$died),]


for (col in 2:length(uniqueds3[1,])){
  uniqueds3[,col] <- as.numeric(as.vector(uniqueds3[,col]))
}

state <- uniqueds2b$state
seas <- uniqueds2b$season
ag <- uniqueds2b$ag
freq <- uniqueds2b$freq
died <- uniqueds2b$died
result <- uniqueds2b$TestResult

testposprop <- uniqueds3

for (row in seq_along(uniqueds3$state)){
  ds3row <- uniqueds3[row,]
  selind1 <- which(state==ds3row$state & seas==ds3row$season & ag==ds3row$ag & died==ds3row$died)
  
  if (length(selind1) > 1){
    pos <- which(result[selind1]==1)
    frq <- uniqueds2b$freq[selind1]
    proppos <- frq[pos]/(frq[pos] + frq[3 - pos])
    testposprop[row,5] <- sum(frq)
    testposprop[row,6] <- as.numeric(proppos)
  } else if (length(selind1)==1 & result[selind1]==0){
    testposprop[row,5] <- freq[selind1]
    testposprop[row,6] <- 0
  } else if (length(selind1)==1 & result[selind1]==1){
    testposprop[row,5] <- freq[selind1]
    testposprop[row,6] <- 1
  }
}

colnames(testposprop) <- c('state','season','ag','died','freq','ppos')
testposprop <- testposprop[which(!is.na(testposprop$freq)),]

testposprop <- testposprop[order(testposprop$state,testposprop$season,testposprop$ag,testposprop$died),]

testposprop$ppos[1:10]
