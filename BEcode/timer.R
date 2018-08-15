tme <- proc.time()
for (i in 1:100000){
  r <- choose(6,2)
}
proc.time() - tme

tme <- proc.time()
for (i in 1:100000){
  r <- bincoeff(6,2)
}
proc.time() - tme