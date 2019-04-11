parms <- c(10,100)
dist <- "Poi"
dist <- "binom"
absc <- c(.05,.15,.8,.9)

abscPrep(absc,parms,dist)


sls <- samples(absc,parms,dist,10000,.75)

hist(sls,100,xlim = c(0,1),freq = F)
lines(ls,yls*100,col = 'red')
