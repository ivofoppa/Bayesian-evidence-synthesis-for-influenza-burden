xls <- 0:10
p0 <- .2
n <- 10
pmeanls <- c(0.08053,0.1656,0.2482,0.3365,0.4162,0.5017,0.5804,0.6686,0.7478,0.8344,0.9163)

xdistls <- sapply(xls,function(x) dbinom(x,n,p0))

sum(pmeanls*xdistls)
