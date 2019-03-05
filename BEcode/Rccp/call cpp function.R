library(Rcpp)

bfolder <- 'C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/'

setwd(paste0(bfolder,'BEcode/Rccp'))

sourceCpp('binom_slice.cpp')
