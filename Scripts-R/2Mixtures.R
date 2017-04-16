setwd('/home/piss/Documents/Extreme/R resources/IRM')
source("1Code_local.R")
load("data1.Rdata")
library(PissoortThesis)

library(evmix)
library(ggplot2)

#hillplot(max_all,try.thresh = c(32,34,36,38,40))
# xi is not >0 here.

fitmweibullgpd(max_all)

mrlplot(max_all)
tcplot(max_all)


hist(above_32$TX)
hist(max_all)
plot(density(max_all))
ggplot(TXTN_closed) + stat_density(aes(x = TX), geom = "line")


dkdengpd(seq(), max_all)




































































