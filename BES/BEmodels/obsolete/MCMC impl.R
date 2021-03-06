library(Rcpp)

model {
  for (seas in 1:5){
    fluhospls[seas] <-  dpois(rfluhospls[seas]*FSNpopls[seas]) ## The # of flu hosp in FSN with not fatal outcome; unobserved
    FSNfluhospls[seas] <-  dbin(ptarr[seas,1],fluhospls[seas]) ## Observed flu hosp. without fatal outcome
    
    rfludeathishls[seas] <- rfludeathls[seas]*(1-poshls[seas])
    fludeathls[seas] <-  dpois(rfludeathishls[seas]*FSNpopls[seas]) ## The total # of flu deaths in FSN; unobserved
    
    FSNfludeathls[seas] <-  dbin(ptarr[seas,2],fludeathls[seas] + (1 - step(fludeathls[seas] - 1))) ## Observed deaths among hospitalized influenza patents
    
    rfluhospls[seas] <-  dunif(0,10)
    rfludeathls[seas] <-  dunif(0,10)
    
    for (k in 1:2){
      ## The three sensitivities, test probs and flu risks correspond to the 3 diff. tests;
      ## the 4th flu risk and 'test prob' represents those not tested
      ptarr[seas,k] <- (pfluarr[seas,k]*ptestarr[k,1,seas]*sens[k,1] + pfluarr[seas,k]*ptestarr[k,2,seas]*sens[k,2] + pfluarr[seas,k]*ptestarr[k,3,seas]*sens[k,3])/(pfluarr[seas,k]*ptestarr[k,1,seas] + pfluarr[seas,k]*ptestarr[k,2,seas] + pfluarr[seas,k]*ptestarr[k,3,seas] + pfluarr[seas,k]*ptestarr[k,4,seas])
      
      nttypearr[k,,seas] <-  dmulti(ptestarr[k,1:4,seas],ntotarr[k,seas])   ### Observed testing freqs
      
      ## priors for testing 
      ptest1[seas,k] <-  dunif(0,1)
      ptestarr[k,1,seas] <- ptest1[seas,k]
      ptest20[seas,k] <-  dunif(0,1)
      ptestarr[k,2,seas] <- ptest20[seas,k]*(1-ptestarr[k,1,seas])
      ptest30[seas,k] <-  dunif(0,1)
      ptestarr[k,3,seas] <- ptest30[seas,k]*(1-ptestarr[k,1,seas]-ptestarr[k,2,seas])
      ptestarr[k,4,seas] <- 1 - ptestarr[k,1,seas] - ptestarr[k,2,seas] - ptestarr[k,3,seas]
      
      pfluarr[seas,k] <-  dunif(0,1)
      
      for (m in 1:3) {
        fluposarr[k,m,seas] <-  dbin(pfluarr[seas,k],nttypearr[k,m,seas] + (1-step(nttypearr[k,m,seas] - 1))) ## actual flu pos/test type, unobserved
        testposarr[k,m,seas] <-  dbin(sens[k,m],fluposarr[k,m,seas]) ## observed flu pos/test type
      }
    }
  }
  ## generating a sensitvity from the given distribution
  for (k in 1:2){
    sens1[k] <-  dnorm(pcrsens[1],1/(pcrsens[2]*pcrsens[2]))
    logsens2[k] <-  dnorm(lrapidsens[1],1/(lrapidsens[2]*lrapidsens[2]))I(,0)
    sens2[k] <- exp(logsens2[k])
    sens3[k] <-  dunif(0,1)
    
    sens[k,1] <- sens1[k]
    sens[k,2] <- sens2[k]
    sens[k,3] <- sens3[k]
  }
}
