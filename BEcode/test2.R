state <- FluSurvdata$state
seas <- FluSurvdata$seas
ag <- FluSurvdata$ag
died <- FluSurvdata$died
TestedFlu <- FluSurvdata$TestedFlu
TestResult <- FluSurvdata$TestResult
TestType <- FluSurvdata$TestType

state <- dataset2$State
seas <- dataset2$season
ag <- dataset2$ag
died <- dataset2$Died
TestedFlu <- dataset2$TestedFlu
TestResult <- dataset2$TestResult
TestType <- dataset2$TestType

selind0 <- which(state=='CA' & seas==1 & ag==5 & died==1 & TestedFlu==1 & TestType==1 & TestResult==0)
selind1 <- which(state=='CA' & seas==1 & ag==5 & died==1 & TestedFlu==1 & TestType==1 & TestResult==1)

length(selind1)/(length(selind0) + length(selind1))

length(which(state=='CA' & seas==1 & ag==5 & died==1  & TestedFlu==1 & TestType==1 & TestResult==0))

sl <- which(state=='CA' & seas==3 & ag==5 & died==1  & TestedFlu==1 & TestType==1 & TestResult==0)
dataset2[sl,]
