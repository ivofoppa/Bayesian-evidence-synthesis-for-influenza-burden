ls <- unique(dataset[which(dataset$Died==1 & (is.na(dataset$ICD9_1) | dataset$ICD9_1!='')),])
length(ls[,1])

ls3 <- sapply(as.vector(ls2),function(x) ifelse(is.character(x) & !is.na(x),substr(x,1,3),'0'),USE.NAMES = F)

ls4 <- sort(unique(ls3))



