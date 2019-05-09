c1 = fupperhull(zval + 1e-6,df,parms,parmind,data,zval,maxind,zind) ;
c2 = fupperhull(zval - 1e-6,df,parms,parmind,data,zval,maxind,zind) ;
c = max(c1,c2) ;

rtio <- 2

(exp(c - mx)/exp(lik(zval,parms,parmind,data) - mx)) > rtio

sapply(1:6, function(s) sum(FSNdata$freq[which(FSNdata$season==s)]))

sapply(1:6, function(s) sum(FSNdata$freq[which(FSNdata$season==s)]))

datasetst <- dataset[which(dataset$State=="CA"),]

datasetcum <- NULL

statels <- as.vector(unique(dataset$State))

for (st in statels) {
  for (seas in 1:6) {
    for (ag in 1:5) {
      ## died
      selseasagind1 <- which(dataset$State==st,dataset$agecat==ag & 
                               dataset$season==seas & dataset$Died==1)
      datasetagseas1 <- dataset[selseasagind1,]
      if (length(selseasagind1) > 0) {
        for (tt in 1:3) {
          numtest1 <- length(which(datasetagseas1$TestedFlu==1 &
                                     datasetagseas1$TestTyperec==tt & datasetagseas1$TestResult==1))
          numtest0 <- length(which(datasetagseas1$TestedFlu==1 &
                                     datasetagseas1$TestTyperec==tt & datasetagseas1$TestResult!=1))
          
          el1 <- c(st,seas,ag,1,1,tt,1,numtest1)
          el0 <- c(st,seas,ag,1,1,tt,0,numtest0)
          datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
        }
      } else {
        el1 <- c(st,seas,ag,1,1,tt,1,0)
        el0 <- c(st,seas,ag,1,1,tt,0,0)
        datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
      }
      num <- length(which(datasetagseas1$TestedFlu!=1))
      el <- c(st,seas,ag,1,0,0,0,num)
      datasetcum <- rbind(datasetcum,el,deparse.level = 0)
      ##did not die
      selseasagind0 <- which(dataset$agecat==ag & 
                               dataset$season==seas & dataset$Died!=1)
      datasetagseas0 <- dataset[selseasagind0,]
      if (length(selseasagind1) > 0) {
        for (tt in 1:3) {
          numtest1 <- length(which(datasetagseas0$TestedFlu==1 &
                                     datasetagseas0$TestTyperec==tt & datasetagseas0$TestResult==1))
          numtest0 <- length(which(datasetagseas0$TestedFlu==1 & 
                                     datasetagseas0$TestTyperec==tt & datasetagseas0$TestResult!=1))
          
          el1 <- c(st,seas,ag,0,1,tt,1,numtest1)
          el0 <- c(st,seas,ag,0,1,tt,0,numtest0)
          datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
        } 
      } else {
        el1 <- c(st,seas,ag,0,1,tt,1,0)
        el0 <- c(st,seas,ag,0,1,tt,0,0)
        datasetcum <- rbind(datasetcum,el1,el0,deparse.level = 0)
      }
      
      num <- length(which(datasetagseas0$TestedFlu!=1))
      el <- c(st,seas,ag,0,0,0,0,num)
      datasetcum <- rbind(datasetcum,el,deparse.level = 0)
    }
  }
}
colnames(datasetcum) <- c('state','season','agecat','died','TestedFlu','TestType','TestResult','freq')
FSNdata <- data.frame(datasetcum)
