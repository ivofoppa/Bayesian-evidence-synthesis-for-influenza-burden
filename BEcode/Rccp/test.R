ls <- seq(0.0001,.999999999,.01)

uhls <- sapply(ls, function(p) fupperhull(p,abscaug,f))

yls <- sapply(ls,fbin)

lhls <- sapply(ls, function(p) flowerhull(p,abscaug,f))


plot(ls,exp(yls),type = 'l',ylim = c(0,.12))
lines(ls,lhls,col = 'blue')
lines(ls,uhls,col = 'red')


plot(ls,uhls,type = 'l')

plot(abscaug,exp(f))

fprob(abscaug,f,zval)
