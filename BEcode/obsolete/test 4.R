ptot <- NULL

ptest10 = ptest10init
ptest20 = ptest20init
ptest30 = ptest30init
ptest40 = ptest40init

ptestarr <- array(0,dim = c(2,4))
sensarr <- array(0,dim = c(2,3))

ptinit <- NULL
for (k in 1:2) {
  ptot[k] <- ptest10[k] + ptest20[k] + ptest30[k] + ptest40[k]
  
  ptestarr[k,1] <- ptest10[k]/ptot[k]
  ptestarr[k,2] <- ptest20[k]/ptot[k]
  ptestarr[k,3] <- ptest30[k]/ptot[k]
  ptestarr[k,4] <- ptest40[k]/ptot[k]
  
  sensarr[k,1] <- sens1init[k]
  sensarr[k,2] <- exp(logsens2init[k])
  sensarr[k,3] <- sens3init[k]
  
  ptinit[k] <- ptestarr[k,1]*sensarr[k,1] + ptestarr[k,2]*sensarr[k,2] + 
    ptestarr[k,3]*sensarr[k,3]
  
}
pt <- ptinit; fluhosp <- fluhospinit

c(pt[1],fluhosp)
c(pt[2],fludeathinit)
c(poshinit,Npideath)
k <- 2; m <- 2

c(pfluinit[k],nttypearr[k,m])
c(sensarrinit[k,m],fluposarrinit[k,m])
