dx <- 0.0001
ls <- seq(0.0001,.999999999,dx)
# ls <- seq(0.2,.7,dx)

uhls <- sapply(ls, function(p) fupperhull(p,abscaug,f,zval,maxind))/(sum(exp(yls))*dx)
# uhls <- sapply(ls, function(p) fupperhull2(p,abscaug,f,zval,maxind))

yls <- sapply(ls,fbin)

lhls <- sapply(ls, function(p) flowerhull(p,abscaug,f))/(sum(exp(yls))*dx)


plot(ls,exp(yls)/(sum(exp(yls))*dx),type = 'l',xlim=c(0,.6),ylim = c(0,max(uhls)))
# plot(ls,exp(yls),type = 'l',xlim=c(.2,.8),ylim = c(0,.2))

lines(ls,lhls,col = 'blue')
lines(ls,uhls,col = 'red')


hist(probls[-c(1:10000)] ,freq = F,xlim=c(0,.6),breaks = 100)
lines(ls,exp(yls)*102,col = 'red')


p1 <- fliksum(abscaug,f,zval,maxind)
p2 <- sapply(seq_along(abscaug[-1]), function(k) sum(sapply(ls[which(ls>=abscaug[k] & ls<abscaug[k + 1])],
                                                             function(x) dx*fupperhull(x,abscaug,f,zval,maxind))))
p1
p2

nsim <- 100000
sls <- NULL

while (length(sls) < nsim) {
  smp <- fpsample(abscaug,lvec,f,zval,maxind)
    sls <- c(sls,smp)
}

hist(sls,breaks = 100, freq  = F,xlim = c(0,1))

pls <- binomslicep(20,100,.3,.1,10000)

