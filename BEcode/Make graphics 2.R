rm(list = ls())
library(ggplot2)
library(scales)

bfolder <- "C:/Users/VOR1/OneDrive/Work related/Misc work/BES and related/Bayesian-evidence-synthesis-for-influenza-burden/"
bfolder2 <- "C:/Users/VOR1/OneDrive/Work related/Misc work/BES and related/Excess mortality/EMrev/" ## project folder

setwd(paste0(bfolder,"BEwriteup/presentations etc."))

source("prelim graphics 2.R")

### Comparing multiplier vs. BE
###################################################################################################
###  Season plots for all age groups, 2012/13 through 2016/17 #####################################
###################################################################################################
hospest <- burdhosarr[which(burdhosarr$Source%in%c("Mult","BE")),]

p11 <- ggplot(hospest,aes(x=Season, y=Estimate, color=Source, group=Source)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, 
                position=position_dodge(0.05)) +
  labs(title = NULL,y = "Estimate") + ylim(0,NA)

p11 + scale_x_discrete(name="Season",labels=c("1"="2012/13","2"= "2013/14","3"= "2014/15","4"= "2015/16","5" ="2016/17")) + 
  scale_color_discrete(name = "Model") + theme(legend.position="top") 

setwd(paste0(bfolder,"BEwriteup/BEmanuscript/BEGraphics"))
ggsave("Hosp_comp_seaon.pdf")

mortest1 <- burddtharr[which(burddtharr$Source%in%c("Mult","BE")),]

p12 <- ggplot(mortest1,aes(x=Season, y=Estimate, color=Source, group=Source)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, 
                position=position_dodge(0.05)) +
  labs(title = NULL,y = "Estimate")+
  ylim(0,NA)

p12 + scale_x_discrete(name="Season",labels=c("1"="2012/13","2"= "2013/14","3"= "2014/15","4"= "2015/16","5" ="2016/17")) + 
  scale_color_discrete(name = "Model") + theme(legend.position = "top")

setwd(paste0(bfolder,"BEwriteup/BEmanuscript/BEGraphics"))
ggsave("Mort_comp_season.pdf")
###################################################################################################
###  Age group plots for season 2016/17 ###########################################################
###################################################################################################
hospest <- burdaghosarr
p31 <- ggplot(hospest,aes(x=agcat, y=Estimate, color=Source, group=Source)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, 
                position=position_dodge(0.05)) +
  labs(title = NULL,y = "Estimate") + ylim(0,NA)

p31 + scale_x_discrete(name="Age Group",labels=c("1"="0-4","2"= "5-17","3"= "18-49",
                                                 "4"= "50-64","5" ="65+")) + 
  scale_color_discrete(name = "Model") + theme(legend.position = "top")
setwd(paste0(bfolder,"BEwriteup/BEmanuscript/BEGraphics"))
ggsave("Hosp_comp_agcat.pdf")

mortest1 <- burdagdtharr[which(burdagdtharr$Source%in%c("Mult","BE")),]

p32 <- ggplot(mortest1,aes(x=agcat, y=Estimate, color=Source, group=Source)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, 
                position=position_dodge(0.05)) +
  labs(title = NULL,y = "Estimate") + ylim(0,NA)

p32 + scale_x_discrete(name="Age Group",labels=c("1"="0-4","2"= "5-17","3"= "18-49",
                                                 "4"= "50-64","5" ="65+")) + 
  scale_color_discrete(name = "Model") + theme(legend.position = "top")
setwd(paste0(bfolder,"BEwriteup/BEmanuscript/BEGraphics"))
ggsave("Mort_comp_agcat.pdf")
