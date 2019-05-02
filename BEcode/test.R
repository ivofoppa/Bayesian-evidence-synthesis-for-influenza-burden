nsims <- 10000
uhls <- NULL
k <- 0
while (k < nsims) {
  k <- k+ 1
  ### Simulated case numbers
  casesls <- controlsls <- NULL
  for (exp in c(0,1)){
    eval(parse(text = paste0('mu <- mu',exp)))
    eval(parse(text = paste0('nu <- nu',exp)))
    cases <- rpois(1,mu)
    casesls[exp + 1] <- cases
    controls <- rpois(1,nu)
    controlsls[exp + 1] <- controls
  }
  data <- data.frame(exp = c(0,1), cases = casesls, controls = controlsls)
  
  for (parmind in 1:2) {
    totList <- superList[[parmind]]
    totList <- abscGroome(totList,data,parmind,15,2)
    
    df <- as.data.frame(totList['df'])
    absc <- df[,1]
    l <- df[,2]
    
    mx <- max(l)
    zval <- unlist(totList['zval'])
    maxind <- unlist(totList['maxind'])
    zind <- unlist(totList['zind'])
    lvec <- fliksum(df,parmind,parms,data,zval,maxind,zind)
    
    superList[[parmind]] <- totList
    psmp <- fpsample(df,parmind,parms,data,zval,maxind,zind,lvec)
    lh <- exp(flowerhull(psmp,df) - mx)
    uh <- exp(fupperhull(psmp,df, parmind, parms, data, zval, maxind, zind) - mx)
    absc <- sort(unique(c(psmp,absc)))
    selind <- which(absc==psmp)
    l1 <- l[1:(selind - 1)]
    l2 <- l[(selind):length(l)]
    l <- c(l1,lik(psmp,parmind,parms,data),l2)
    
    df <- data.frame(absc,l)
    totList[["df"]] <- df
    
    totList <- abscGroome(totList,data,parmind,15,2)
    
    maxind <- maxindC(df, parmind, parms, data) ;
    zval <- intsct(df,parmind,parms,data,maxind) ;
    zind <- fzind(df,parmind,parms,data,zval) ;
    
    lvec <- fliksum(df,parmind,parms,data,zval,maxind,zind)
    
    totList[["df"]] <- df
    superList[[parmind]] <- totList
    }
}
    acc <- 0
    while (acc==0) {
      uran <- runif(1)
      psmp <- fpsample(df,parmind,parms,data,zval,maxind,zind,lvec)
      
      lh <- exp(flowerhull(psmp,df) - mx)
      uh <- exp(fupperhull(psmp,df, parmind, parms, data, zval, maxind, zind) - mx)
      
      if (uran < lh/uh) {
        acc <- 1
      } else {
        lval <- exp(lik(psmp,parmind,parms,data) - mx)
        if (uran < lval/uh) {
          acc <- 1
        } else if (lval/uh < crit) {
          absc <- sort(unique(c(psmp,absc)))
          selind <- which(absc==psmp)
          l1 <- l[1:(selind - 1)]
          l2 <- l[(selind):length(l)]
          l <- c(l1,lik(psmp,parmind,parms,data),l2)
          
          df <- data.frame(absc,l)
          totList[["df"]] <- df
          
          totList <- abscGroome(totList,data,parmind,15,2)
          
          maxind <- maxindC(df, parmind, parms, data) ;
          zval <- intsct(df,parmind,parms,data,maxind) ;
          zind <- fzind(df,parmind,parms,data,zval) ;
          
          lvec <- fliksum(df,parmind,parms,data,zval,maxind,zind)
          
          totList[["df"]] <- df
          superList[[parmind]] <- totList
        } 
      }
    }
    parms[parmind] <- psmp
    superList[[1]][["parms"]] <- parms
    superList[[2]][["parms"]] <- parms
  }
  mcmcarr <- rbind(mcmcarr,parms, deparse.level = 0)
  if (length(mcmcarr[,1])%%1000==0) {
    cat(length(mcmcarr[,1]),' iterations!\n')
  }
}

