burdagdtharr <- array(0,dim = c(5*5,5))
burdaghosarr <- array(0,dim = c(5*2,5))

for (ag in 1:5){
  selind <- which(burdendata0$X==ag & burdendata0$season=="2016-2017")
  burdagdtharr[ag,] <- c("Mult", ag,as.numeric(burdendata0$dcount[selind]),as.numeric(burdendata0$dlower[selind]),as.numeric(burdendata0$dupper[selind]))
  burdaghosarr[ag,] <- c("Mult", ag,as.numeric(burdendata0$hcount[selind]),as.numeric(burdendata0$hlower[selind]),as.numeric(burdendata0$hupper[selind]))
}

setwd(paste0(bfolder,'BEmcmc'))

infile <- "codafileFSN_OSH2012-17.RData"
load(infile)

for (ag in 1:5){
  codals <- codaagList[[ag]]
  BEdth <- round(as.numeric(quantile(codals[,5],probs = c(.5,.025,.975))))
  BEhos <- round(as.numeric(quantile(codals[,5 + 5],probs = c(.5,.025,.975))))
  burdagdtharr[ag + 5,] <- c("BE",ag,BEdth)
  burdaghosarr[ag + 5,] <- c("BE",ag,BEhos)
}

burdaghosarr <- data.frame(burdaghosarr)
colnames(burdaghosarr) <- c("Source","agcat","Estimate","lower","upper")

for (col in 2:length(burdaghosarr[1,])) {
  burdaghosarr[,col] <- as.numeric(as.vector(burdaghosarr[,col]))
}
burdaghosarr$agcat <- as.factor(burdaghosarr$agcat)
###################################################################################################
###  Reading-in EM estimates ######################################################################
###################################################################################################
fromyr <- 2012
toyr <- 2017
fromwk <- 40
towk <- 23
nseas <- toyr-fromyr
niter <- 1000

setwd(paste0(bfolder2,"EMrevoutput"))
qualls <- c("flu","vir","vir alt")
for (ind in 1:3) {
  eval(parse(text = paste0("fname <- paste0(\"codaarr ",qualls[ind], " \",",niter,",\" National \",",fromyr,",\"-\",",substr(toyr,3,4),",\" RC.RData\")")))
  load(fname)
  EMcoda <- codaarrag
  for (ag in 1:5) {
    EMcodasel <- EMcoda[[5]]
    EMest <- round(as.numeric(quantile(EMcodasel[,5],probs = c(.5,.025,.975))))
    burdagdtharr[(ind + 1)*5 + ag,1] <- qualls[ind]
    burdagdtharr[(ind + 1)*5 + ag,2:5] <- as.numeric(c(ag,EMest))
  }
}

burdagdtharr <- data.frame(burdagdtharr)
colnames(burdagdtharr) <- c("Source","agcat","Estimate","lower","upper")

for (col in 2:length(burdagdtharr[1,])) {
  burdagdtharr[,col] <- as.numeric(as.vector(burdagdtharr[,col]))
}

burdagdtharr$agcat <- as.factor(burdagdtharr$agcat)
setwd(paste0(bfolder,"BEwriteup"))
###################################################################################################
