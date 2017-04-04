library(evd)
library(fExtremes)
library(ismev)
library(evir)
library(RColorBrewer)
library(AoE)

# TX_GEV <- gev(TXTN_IRM_closed$TX,block=365)
# TN_GEV <- gev(TXTN_IRM_closed$TN,block = 365)
# plot(TX_GEV)


##############################################
#############################"################
# Remember to faire des espace entre chaque , = + etc etc pour faiciliter la visibilit?

#==============================================================
############ First try : standard plots ######################
#-------------------------------------------------------------

x=seq(-5,5,length=5000)

mu <- 0.5
sig <- 0.5
#################    Attention frechet scaled pasbonne !! 
xi_gumbel <- 0
gumbel_dens <-  exp(-x) * exp(-exp(-x))

xi_frechet05 <- 0.5
s <- (1 + xi_frechet05 * (x - mu)/sig)^(-(xi_frechet05)^-1 - 1)
t <- (1 + xi_frechet05 * (x - mu)/sig)^(-(xi_frechet05)^-1)
u <-  (x - mu)/sig >  -1/xi_frechet05
frechet_dens05 <- sig^-1 *  s * exp(-t) * u
xi_weibull05 <- -0.5
s <- (1 + xi_weibull05 * (x - mu)/sig)^(-(xi_weibull05)^-1 -1)
t <- (1 + xi_weibull05 * (x - mu)/sig)^(-(xi_weibull05)^-1)
u <- (x - mu)/sig  < -1/xi_weibull05
weibull_dens05 <- s * exp(-t) * u

xi_frechet015 <- 0.15
s <- (1 + xi_frechet015 * (x - mu)/sig)^(-(xi_frechet015)^-1 - 1)
t <- (1 + xi_frechet015 * (x - mu)/sig)^(-(xi_frechet015)^-1)
u <-  (x - mu)/sig >  -1/xi_frechet015
frechet_dens015 <- sig^-1 *  s * exp(-t) * u
xi_weibull015 <- -0.15
s <- (1 + xi_weibull015 * (x - mu)/sig)^(-(xi_weibull015)^-1 -1)
t <- (1 + xi_weibull015 * (x - mu)/sig)^(-(xi_weibull015)^-1)
u <- (x - mu)/sig  < -1/xi_weibull015
weibull_dens015 <- s * exp(-t) * u

xi_frechet1 <- 0.6
s <- (1 + xi_frechet1 * (x - mu)/sig)^(-(xi_frechet1)^-1 - 1)
t <- (1 + xi_frechet1 * (x - mu)/sig)^(-(xi_frechet1)^-1)
u <-  (x - mu)/sig >  -1/xi_frechet1
frechet_dens1 <- sig^-1 *  s * exp(-t) * u
xi_weibull1 <- -0.6
s <- (1 + xi_weibull1 * (x - mu)/sig)^(-(xi_weibull1)^-1 -1)
t <- (1 + xi_weibull1 * (x - mu)/sig)^(-(xi_weibull1)^-1)
u <- (x - mu)/sig  < -1/xi_weibull1
weibull_dens1 <- s * exp(-t) * u

###################
plot(x, frechet_dens015, type="l", col="blue", ylim=c(0,0.5),main = "Generalized extreme value densities",
     ylab = "Density", xlab = "x",lwd=2,sub = "all with mu=0, sigma=1")
lines(x, weibull_dens05, col="green", lwd=2)
lines(x,gumbel_dens, col="red", lwd=2)
legend(-4,0.5,c("xi = -0.5","xi=0","xi=+0.5"),lwd=2,col=c("green","red","blue"))
#######################



#==============================================================
#################### ggplots #################################
#-------------------------------------------------------------


library(ggplot2)
library(gganimate)   #!!!!!

GEV_g <- as.data.frame(cbind(x,gumbel_dens,xi_gumbel))
GEV_g$xi_gumbel <- as.factor(GEV_g$xi_gumbel)

GEV_f015 <- as.data.frame(cbind(x,frechet_dens015,xi_frechet015))
GEV_w015 <- as.data.frame(cbind(x,weibull_dens015,xi_weibull015))

GEV_f05 <- as.data.frame(cbind(x,frechet_dens05,xi_frechet05))
GEV_w05 <- as.data.frame(cbind(x,weibull_dens05,xi_weibull05))

GEV_f1 <- as.data.frame(cbind(x,frechet_dens1,xi_frechet1))
GEV_w1 <- as.data.frame(cbind(x,weibull_dens1,xi_weibull1))


colnames(GEV_g) <- c("x",'Density',expression(xi))
colnames(GEV_f015) <- colnames(GEV_w015) <- c("x",'Density',expression(xi))
colnames(GEV_f02) <- colnames(GEV_w02) <- c("x",'Density',expression(xi))
colnames(GEV_f05) <- colnames(GEV_w05) <- c("x",'Density',expression(xi))
colnames(GEV_f1) <- colnames(GEV_w1) <- c("x",'Density',expression(xi))

GEV_densities <- rbind(GEV_g,GEV_w015,GEV_f015,#GEV_f05,GEV_w05,
                       GEV_f1,GEV_w1)

#GEV_densities <- as.data.frame(GEV_densities)

GEV_densities$xi <- as.factor(GEV_densities$xi) 

####   try
# ggplot(GEV_densities,aes(x=x,y=density,fill=xi)) +
 # geom_line(aes(fill=xi), size=2,show.legend = TRUE) +
  # scale_colour_gradient(low="red") + pres + theme_minimal()


pres <- labs(title=expression(paste('Generalized Extreme Value density (?=0,',sigma,"=1)")),
             colour=expression(paste(xi,"="),linetype=expression(paste(xi,"="))))
pres1 <- labs(title=expression(paste("With decreased shape parameters ",xi)),
             colour=expression(paste(xi,"=")))
pres2 <- labs(title=expression(paste("With increased shape parameters ",xi)),
             colour=expression(paste(xi,"=")))

plot1 <- rbind(GEV_f015,GEV_f05) 

plot2 <- rbind(GEV_w02,GEV_w1) 
plot2$xi <- as.factor(plot2$xi)

plot05 <- rbind(GEV_w05,GEV_w1) 
plot05$xi <- as.factor(plot2$xi)

plot2 <- rbind(GEV_w015,GEV_w05) 
plot2$xi <- as.factor(plot2$xi)

GEV_g <- as.data.frame(cbind(x,dgumbel(x),rep(0,length(x))))
GEVf015 <-as.data.frame(cbind(x,dfrechet(x,shape=0.15),rep(0.15,length(x))))
GEVf05 <- as.data.frame(cbind(x,dfrechet(x,shape=0.5),rep(0.5,length(x))))

colnames(GEVf015) <- colnames(GEVf05) <- colnames(GEV_g)<- c("x",'Density',expression(xi))

GEVw015 <- 
GEVw05 <- 

GEV_densities <- GEV_g
GEV_densities <- cbind.data.frame(GEV_g,GEV_w015,GEVf015,GEVf05,GEV_w05)
                         #GEVf1,GEV_w1)
GEV_densities$xi <- as.factor(GEV_densities$xi) 
colnames(GEV_densities) <- c("x",'Density',"xi1",............)
GEV_densities <- as.data.frame(GEV_densities)

# http://stackoverflow.com/questions/14604435/turning-off-some-legends-in-a-ggplot
# http://stackoverflow.com/questions/20378276/legend-does-not-show-line-type-in-ggplot2-density-plot
# http://stackoverflow.com/questions/23635662/editing-legend-text-labels-in-ggplot
ggplot(data=GEV_densities) + 
  geom_line(aes(x=x,y=GEV_densities[,2],linetype="0",color="0"),size=1.1) +
  geom_line(aes(x=x,y=GEV_densities$frechet_dens015,linetype="0.15",color="0.15"),size=1.1) +
  geom_line(aes(x=x,y=GEV_densities[,11],linetype="0.5",color="0.5"),size=1.1)  +
  geom_line(aes(x=x,y=GEV_densities[,5],linetype="-0.15",color="-0.15"),size=1.1) + 
  geom_line(aes(x=x,y=GEV_densities[,14],linetype="-0.5",color="-0.5"),size=1.1) +
  scale_linetype_manual(values=c("0"=1, "0.15"=4, "0.5"=1,
                                 "-0.15"=4, "-0.5"=1), name="xi") +
  scale_colour_manual(values=c("0"="black","0.15"="blue", "0.5"="blue",
                                "-0.15"="red", "-0.5"="red"), name="xi") +
  theme_bw() + labs(title=expression(paste('Generalized Extreme Value density (?=0,',sigma,"=1)")))



ggplot(data=GEV_densities,fill=xi) + 
  geom_line(data=GEV_g,aes(x=x,y=Density,linetype="0",color="0"),size=1.1) +
  geom_line(data=GEV_f01,aes(x=x,y=Density,linetype="0.15",color="0.15"),size=1.1) +
  geom_line(data=GEV_f05,aes(x=x,y=Density,linetype="0.5",color="0.5"),size=1.1)  +
  geom_line(data=GEV_w015,aes(x=x,y=Density,linetype="-0.15",color="-0.15"),size=1.1) + 
  geom_line(data=GEV_w05,aes(x=x,y=Density,linetype="-0.5",color="-0.5"),size=1.1) +
  scale_linetype_manual(values=c("0"=1, "0.15"=4, "0.5"=1,
                                 "-0.15"=4, "-0.5"=1), name="xi") +
  scale_colour_manual(values=c("0"="black","0.15"="blue", "0.5"="blue",
                               "-0.15"="red", "-0.5"="red"), name="xi") +
  theme_bw() + labs(title=expression(paste('Generalized Extreme Value density (?=0,',sigma,"=1)")))



ggplot(data=GEV_densities,col=xi,linetype=xi) + geom_line(data=GEV_g,aes(x=x,y=Density),size=1.1) +
  geom_line(data=plot1,aes(x=x,y=Density,linetype=xi,color="darkred"),size=1.1) +
  geom_line(data=plot2,aes(x=x,y=Density,linetype=xi,color="darkblue"),size=1.1)  +
  scale_linetype_manual(values=c("0.5"="dotted","1"="solid","-0.2"="dotted","-1"="solid"),
                        name="sha") +
  scale_colour_manual("0.15"="blue",,name=expression(xi)) +theme_bw() +
  labs(title=expression(paste('Generalized Extreme Value density (?=0,',sigma,"=1)")))
       
       
  ggplot(data=GEV_densities,col=xi,linetype=xi) + geom_line(data=GEV_g,aes(x=x,y=Density),size=1.1) +
    geom_line(data=plot1,aes(x=x,y=Density,linetype=xi,color="darkred"),size=1.1) +
    geom_line(data=plot2,aes(x=x,y=Density,linetype=xi,color="darkblue"),size=1.1)  +
    scale_linetype_manual(values=c("solid","longdash","solid","longdash"),
                          guide=guide_legend(title=expression(xi))) +theme_bw() +
    
    labs(title=expression(paste('Generalized Extreme Value density (?=0,',sigma,"=1)")))
  

ggplot(data = GEV_densities,aes(x=x,y=Density,fill=xi,color=xi,linetype=xi)) + 
  geom_line(lwd=1,size=1.1) +
  scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 2))) +
  scale_color_manual(values = c(brewer.pal(4, "Set3"), brewer.pal(2, "Set3"))) +
  opts(title = "BC mortality") +
   + theme_minimal() 



ggplot() + 
  geom_line(data = GEV_densities,aes(x=x,y=Density,fill=xi,color=xi),size=1.1) +
  coord_cartesian(ylim = c(0, 0.57)) + pres1 + theme_minimal() 
ggplot() + 
  geom_line(data = GEV_densities,aes(x=x,y=Density,fill=xi,color=xi),size=1.1) +
  pres2 + theme_minimal() 


gumbel_auto <- dgumbel(x,loc =1)
frechet_auto <- dfrechet(x,loc=1)
weibull_auto <- drweibull(x,loc=1)

gev1_auto <-cbind(x,gumbel_auto,xi_gumbel)
gev2_auto <-  cbind(x,frechet_auto,xi_frechet02)
     gev3_auto <-    cbind(x,weibull_auto,xi_weibull02)
     
GEV_dens_auto <- as.data.frame(rbind(gev1_auto,gev2_auto,gev3_auto))
colnames(GEV_dens_auto) <- c("x",'Density',expression(xi))
GEV_dens_auto$xi <- as.factor(GEV_dens_auto$xi) 

pres <- labs(title=expression(paste('Generalized Extreme Value density (?=0,',sigma,"=1)")),
             colour=expression(paste(xi,"=")))

ggplot() + 
  geom_line(data = GEV_dens_auto,aes(x=x,y=Density,fill=xi,color=xi,linetype=xi),size=1.1) +
  scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 2))) + 
  pres + theme_minimal() 


# as we will mainly use Fr?chet dist ? we represent the Frechet density 
#   for various values of parameters

frechet_auto <- dfrechet(x,shape=1.5)
frechet_auto <- dfrechet(x,loc=1)
frechet_auto <- dfrechet(x,loc=1)


ggplot() + 
  geom_line(data = GEV_dens_auto,aes(x=x,y=Density,fill=xi,color=xi),size=1.1) +
  pres + theme_minimal() 


# http://tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.html milieu bas page



library(POT)



