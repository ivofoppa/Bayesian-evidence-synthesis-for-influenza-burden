dx <- 0.0001
ls <- seq(0.0001,.999999999,dx)

uhls <- sapply(ls, function(p) fupperhull(p,abscaug,f,zval))

yls <- sapply(ls,fbin)

lhls <- sapply(ls, function(p) flowerhull(p,abscaug,f))


plot(ls,exp(yls),type = 'l',ylim = c(0,.15))
lines(ls,lhls,col = 'blue')
lines(ls,uhls,col = 'red')


p1 <- fprob(abscaug,f,zval)
p2 <- sapply(seq_along(abscaug[-1]), function(k) sum(sapply(ls[which(ls>=abscaug[k] & ls<abscaug[k + 1])],
                                                             function(x) dx*fupperhull2(x,abscaug,f,zval))))

p1
p2

fupperhull2(zval+.00002,abscaug,f,zval)

pls <- seq(0,.1,.001) - .05 + zval

sapply(pls, function(p) fupperhull2(p,abscaug,f,zval))

ls <- NULL 
for (k in 1:100) {
  ls[k] <- fpsample(abscaug,lvec,f,zval)
}
mean(ls)
