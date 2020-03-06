ag <- 5

sapply(1:nseas,function(s) sum(OSHcumdata$pi[which(OSHcumdata$agecat==ag & OSHcumdata$season==s & OSHcumdata$osh==1)])/
         sum(OSHcumdata$pi[which(OSHcumdata$agecat==ag & OSHcumdata$season==s)]))

sapply(1:nseas,function(s) sum(OSHcumdata$rcu[which(OSHcumdata$agecat==ag & OSHcumdata$season==s & OSHcumdata$osh==1)])/
         sum(OSHcumdata$rcu[which(OSHcumdata$agecat==ag & OSHcumdata$season==s)]))
