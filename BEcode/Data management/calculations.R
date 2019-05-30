sapply(1:8,function(s) 1 - nttypearr[1,4,s]/sum(nttypearr[1,,s]))
sapply(1:8,function(s) 1 - nttypearr[2,4,s]/sum(nttypearr[2,,s]))

sapply(1:8,function(s) 1 - sum(nttypearr[,4,s])/sum(nttypearr[,,s]))

selindst <- which(dataset$TestedFlu==1)
dataset2 <- dataset[selindst,]
unique(dataset$State)
