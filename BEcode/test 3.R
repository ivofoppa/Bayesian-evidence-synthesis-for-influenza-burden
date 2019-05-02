pt <- c(.5,.2,.1,.2)
nttype <- c(961,267,37)
sens <- c(.9,.7,.6,0)
testpos <- c(264,65,7)
pos <- round(testpos/sens)
sum(testpos)/sum(pos)
pflu <- 0.84

N <- 10000 ## No. of all flu in FSN
NFSN <- 6000 ## No. of detected flu

P1 <- 6300 ## No. positive tested with test 1
TP1 <- 5670

## Lik(pt1|...)
pt1ls <- seq(0,1 - pt2,.0001)
likpt1ls <- sapply(pt1ls, function(pt1) dbinom(NFSN,N,pt%*%sens))
plot(pt1ls,likpt1ls,type = 'l')

## Lik(sens1|...)
sens1ls <- seq(0,1,.0001)
liksens1ls <- sapply(sens1ls, function(sens1) dbinom(NFSN,N,pt%*%replace(sens,1,sens1)) * dbinom(TP1,P1,sens1))
plot(sens1ls,liksens1ls,type = 'l')

## Lik(pflu|...)
pfluls <- seq(.1,.6,.001)
likpfluls <- sapply(pfluls, function(pflu) dbinom(testpos[1],nttype[1],sens[1]*pflu) *
                      dbinom(testpos[2],nttype[2],sens[2]*pflu) *
                      dbinom(testpos[3],nttype[3],sens[3]*pflu))

