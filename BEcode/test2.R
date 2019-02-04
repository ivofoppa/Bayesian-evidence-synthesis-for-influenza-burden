state <- FluSurvdata$state
seas <- FluSurvdata$seas
ag <- FluSurvdata$ag
died <- FluSurvdata$died
TestedFlu <- FluSurvdata$TestedFlu
TestResult <- FluSurvdata$TestResult
TestType <- FluSurvdata$TestType

state <- dataset3$State
seas <- dataset3$season
ag <- dataset3$ag
died <- dataset3$Died
TestedFlu <- dataset3$TestedFlu
TestResult <- dataset3$TestResult
TestType <- dataset3$TestType

selind0 <- which(state=='NY' & seas==1 & ag==5 & died==1 & TestedFlu==1 & TestType==0 & TestResult==0)
selind1 <- which(state=='NY' & seas==1 & ag==5 & died==1 & TestedFlu==1 & TestType==0 & TestResult==1)

state <- as.vector(FluSurvdatashort$state)
seas <- as.numeric(as.vector(FluSurvdatashort$seas))
ag <- as.numeric(as.vector(FluSurvdatashort$ag))
died <- as.numeric(as.vector(FluSurvdatashort$died))
TestedFlu <- as.numeric(as.vector(FluSurvdatashort$TestedFlu))
TestResult <- as.numeric(as.vector(FluSurvdatashort$TestResult))


selind0 <- which(state=='CA' & seas==3 & ag==5 & died==1 & TestedFlu==1 & TestResult==0)
selind1 <- which(state=='CA' & seas==3 & ag==5 & died==1 & TestedFlu==1 & TestResult==1)

FluSurvdatashort$freq[selind0]/(FluSurvdatashort$freq[selind0] + FluSurvdatashort$freq[selind1])
