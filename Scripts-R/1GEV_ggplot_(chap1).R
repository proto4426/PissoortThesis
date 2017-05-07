# *Antoine Pissoort* # 
## Script providing the main plot for the thesis in the chapter 1 : GEV ## 
##########################################################################

library(evd)
library(fExtremes)
library(ismev)
library('ggplot2')


## Create data frames for ggplot 

'GEV.df.Fun' <-
  function (x = seq(-10, 10, length = 10^3), mu1 = 0, sig = 1, ksi = 0) {
   df <- data.frame(x = x, density = evd::dgev(x, loc = mu1, scale = sig, shape = ksi),
                    xi = as.factor(ksi), mu = as.factor(mu1), scale = as.factor(sig))
 return(df)
}

'GEVdfFun' <- 
  function (x = seq(-10, 10, length = 10^4), mu = 0, sig = 1, ksi = 0) {
    if (ksi ==0) { dens <-  exp(-x) * exp(-exp(-x)) }
    else 
    s <- (1 + ksi * (x - mu)/sig)^(-(ksi)^-1 - 1)
    t <- (1 + ksi * (x - mu)/sig)^(-(ksi)^-1)
    if (ksi < 0) {dens <-  s * exp(-t) * ( (x - mu)/sig  < -1/ksi ) }
    if (ksi > 0) {dens <- sig^{-1} * s * exp(-t) * ( (x - mu)/sig  > -1/ksi ) }
    
 df <- data.frame(x = x, density = dens,xi = as.factor(ksi), 
                  mu = as.factor(mu), scale = as.factor(sig))
return(df)
}



GEVdf <- rbind(GEV.df.Fun(ksi = -.5),GEV.df.Fun(ksi = 0),GEV.df.Fun(ksi = .5))
GEVdf <- rbind(GEVdfFun(ksi = -.5),GEVdfFun(ksi = 0),GEVdfFun(ksi = .5))

endpoint_w <- 0 - (1 / -0.5)
endpoint_f <- 0 - (1 / 0.5)
dens_f <- ifelse(GEVdf[GEVdf$xi == 0.5,]$density < endpoint_f, NA,
                GEVdf[GEVdf$xi == 0.5,]$density )
GEVdf[GEVdf$xi == 0.5,]$density <- dens_f


# plot the normal distribution as reference 

GEVdf <- cbind(GEVdf, norm = dnorm(GEVdf$x))



GEVdf[GEVdf$density < 10^{-312}, ]$density <- NA
pres <- labs(title = expression(paste(underline(bold('Generalized Extreme Value density (')),underline('µ'),
                                    underline(bold('=0,')),
                                    underline(bold(sigma)), underline(bold("=1)")))),
                              colour = expression(paste(xi,"="),linetype=expression(paste(xi,"="))))
pres <- labs(title = expression(paste(underline(bold('Generalized Extreme Value density')), " (µ",
                                    '=0,',
                                    sigma, "=1)")),
             colour = expression(paste(xi,"=")),linetype = expression(paste(xi,"=")))


gf <- ggplot(GEVdf, aes(x = x, y = density, colour = xi )) + 
  geom_line() + pres +  coord_cartesian(xlim = c(-8, 8)) + 
  geom_line(aes(x = x, y = norm, col = "normal"), col = "black", linetype = 3)+
  theme_piss(20, 15, theme_classic() ) +
  theme(legend.title = element_text(colour="#33666C", 
                                    size=18, face="bold")) +
  theme(legend.key = element_rect(colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  geom_point(aes(x = endpoint_f, y = 0),size = 3.5) + 
  geom_point(aes(x = endpoint_w, y = 0), col="red",size = 3.5)
gf

gright <- ggplot(GEVdf, aes(x = x, y = density, colour = xi )) + 
  geom_line() + pres + 
  geom_line(aes(x = x, y = norm, col = "normal"), col = "black", linetype = 3) +
  coord_cartesian(xlim = c(1.8, 5.5), ylim = c(0,.15)) +
  theme_piss(18, 15, theme_minimal()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  theme(axis.line = element_line(color="#33666C", size = .3)) +
  theme(legend.title = element_text(colour="#33666C", 
                                    size=18, face="bold")) + 
  theme(axis.title = element_blank()) +
  labs(title = expression(paste(bold("Right tails")))) +
  geom_point(aes(x = endpoint_w, y = 0), col = "red", size = 3) +
  scale_color_discrete(guide=F) 
  
gright

gleft <- ggplot(GEVdf, aes(x = x, y = density, colour = xi )) +
  geom_line() + pres + 
  geom_line(aes(x = x, y = norm, col = "normal"), col = "black", linetype = 3) +
  coord_cartesian(xlim = c(-0.9,-4.5), ylim = c(0,.2)) +
  theme_piss(18, 15, theme_minimal()) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color="#33666C", size = .3),
        legend.title = element_text(colour="#33666C", 
                                    size = 18, face="bold")) + 
  theme(axis.title = element_blank()) +
  labs(title = expression(paste(bold("Left tails")))) +
  geom_point(aes(x = endpoint_f, y = 0), size = 3) + 
  scale_color_discrete(guide=F) 
gleft

gleft2 <- ggplot(GEVdf, aes(x = x, y = density, colour = xi )) +
  geom_line() + 
  geom_line(aes(x = x, y = norm, col = "normal"), col = "black", linetype = 3) +
  coord_cartesian(xlim = c(-1.2,-4.5), ylim = c(0.0001,.01)) +
  theme_piss(18, 15, theme_minimal()) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color="#33666C", size = .3),
        legend.title = element_text(colour="#33666C", 
                                    size = 18, face="bold")) + 
  theme(axis.title = element_blank()) +
  labs(title = expression(paste("Zoom : ", bold("Left tails")))) +
  geom_point(aes(x = endpoint_f, y = 0), size = 3) + 
  scale_color_discrete(guide=F) 
gleft2


library(grid)
vpleft <- viewport(width = 0.29, 
                height = 0.43, 
                x = 0.234, 
                y = 0.42)
vpright <- viewport(width = 0.28, 
                   height = 0.415, 
                   x = 0.73, 
                   y = 0.43)
vpleft2 <- viewport(width = 0.29, 
                   height = 0.13, 
                   x = 0.234, 
                   y = 0.12)
gf 
print(gleft2, vp = vpleft)
print(gright, vp = vpright)
#print(gleft2, vp = vpleft2)





## Inspect the influence of the parameters mu and sigma on the shape of the distribution


x <- seq(min(TXTN_closed$TX), max(TXTN_closed$TX) +4, length = 10^3)
GEVmu <- rbind(GEV.df.Fun(x = x, ksi = -.5, mu = 25),GEV.df.Fun(x = x,ksi = 0, mu = 25),
               GEV.df.Fun(x = x,ksi = .5, mu = 25),
               GEV.df.Fun(x = x,ksi = -.5, mu = 30),GEV.df.Fun(x = x,ksi = 0, mu = 30),
               GEV.df.Fun(x = x,ksi = .5, mu = 30), 
               GEV.df.Fun(x = x,ksi = -.5, mu = 35),GEV.df.Fun(x = x,ksi = 0, mu = 35),
               GEV.df.Fun(x = x,ksi = .5, mu = 35))
titles <- labs( title = expression(paste("GEV densities where ", sigma, "=1 and different values for ", mu)),
                caption = expression(italic("The endpoints of the distributions are not clearly displayed here for convenience, but you can easily find them.")),
                colour = expression(paste(xi,"=")),linetype = expression(paste(mu,"=")))
gmu <- ggplot(GEVmu, aes(linetype = mu, colour = xi)) + geom_line(aes(x = x, y = density)) + 
  theme_bw() + theme_piss(legend.position = c(.96, .58)) + coord_cartesian(xlim = c(22, 40), ylim = c(0.015, .435)) +
  geom_vline(xintercept = 25) +   geom_vline(xintercept = 30, linetype = 3) + 
  geom_vline(xintercept = 35, linetype = 2) + 
  guides(colour = guide_legend(override.aes = list(size = 1.5))) + titles
gmu

GEVsig <- rbind(GEV.df.Fun(ksi = -.5, sig = .5),GEV.df.Fun(ksi = 0, sig = .5),
                GEV.df.Fun(ksi = .5, sig = .5),
                GEV.df.Fun(ksi = -.5, sig = 1),GEV.df.Fun(ksi = 0, sig = 1),
                GEV.df.Fun(ksi = .5, sig = 1), 
                GEV.df.Fun(ksi = -.5, sig = 1.5),GEV.df.Fun(ksi = 0, sig = 1.5),
                GEV.df.Fun(ksi = .5, sig = 1.5))
# levels(GEVsig$scale) <- c(expression(paste(sigma, "= 0.5")), expression(paste(sigma, "= 1")),
#                           expression(paste(sigma, "= 1.5" )))
titles <- labs( title = expression(paste("GEV densities where  µ=0 and ", sigma, " =")), 
                colour = expression(paste(xi,"=")))
gsig <- ggplot(GEVsig, aes(colour = xi)) +
  geom_vline(xintercept = 0) +
  geom_line(aes(x = x, y = density)) + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(0.03,.81)) + theme_piss(18, 15, legend.position = c(.96, .79)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) + 
  facet_wrap(~scale) + titles
gsig


#grid_arrange_legend(gmu, gsig, ncol = 1, nrow = 2, position = "bottom")
grid.arrange(gmu, gsig, nrow = 2)




#### For our fitted model ? 
gev_data <- GEVdfFun(mu = gev_fit$mle[1], sig = gev_fit$mle[2], 
                     ksi = gev_fit$mle[3] )
gev_data <- GEV.df.Fun(mu1 = gev_fit$mle[1], sig = gev_fit$mle[2], 
                       ksi = gev_fit$mle[3] )
ggplot(gev_data) +
  geom_line(aes(x = x, y = density)) + coord_cartesian(xlim = c(-5, 5)) +
  theme(plot.title = element_text(size=18, hjust=0.5, 
                                  colour = "#330000", face="bold")) +
  theme(axis.title = element_text(face = "bold", size= 15,
                                  colour = "#330000")) +
  theme(legend.title = element_text(colour="#330000", 
                                    size=18, face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +theme_bw()





################# For Beamer Presentation 

p <- ggplot(TXTN_closed) + geom_density(aes(x = TX), fill = "grey") + 
  labs(title = "Distribution of maximum daily temperatures (TX)") 

d <- ggplot_build(p)$data[[1]]

p + geom_area(data = subset(d, x > 30), aes(x=x, y=y), fill="red") +
  geom_segment(x=30, xend = 30,
               y=0, yend=approx(x = d$x, y = d$y, xout = 30)$y,
               colour="blue", size=2) + theme_classic() +
  theme(plot.title = element_text(size=18, hjust=0.5, 
                                  colour = "#33666C", face="bold")) +
  theme(axis.title = element_text(face = "bold", size= 15,
                                  colour = "#33666C")) +
  theme(legend.title = element_text(colour="#33666C", 
                                    size=18, face="bold")) +
  geom_segment(aes(x = 30, xend = 31, y = .01, yend = .0029),
               arrow = arrow(length = unit(0.5, "cm")), col = "red") + 
  geom_segment(aes(x = 25, xend = 29.5, y = .004, yend = .002),
               arrow = arrow(length = unit(0.8, "cm")), col = "blue") +
  annotate("text", label = "Region of \n interest",
           x = 28, y = .015, col = "red", size = 9, fontface = 2) + 
  annotate("text", label = "Threshold(?)",
           x = 18, y = .007, col = "blue", size = 8, fontface = 1) 




library(evd)

'gpd.dfFun' <-
  function (x = seq(0, 100, length = 10^3), mu1 = 0, sig = 1, ksi = 0) {
    df <- data.frame(x = x, density = evd::dgpd(x, loc = mu1, scale = sig, shape = ksi),
                     xi = as.factor(ksi), mu = as.factor(mu1), scale = as.factor(sig))
    return(df)
  }

GPDdf <- rbind(gpd.dfFun(ksi = -.2),gpd.dfFun(ksi = 0),gpd.dfFun(ksi = 1))



gf <- ggplot(GPDdf, aes(x = x, y = density, colour = xi )) + 
  geom_line() +
  pres +  coord_cartesian(xlim = c(0.5, 6)) + theme_classic() +
  #geom_line(aes(x = x, y = norm), linetype = 3)+
  theme(plot.title = element_text(size = 20, hjust=0.5, 
                                  colour = "#33666C", face="bold")) +
  theme(axis.title = element_text(face = "bold", size= 15,
                                  colour = "#33666C")) +
  theme(legend.title = element_text(colour="#33666C", 
                                    size=18, face="bold")) + 
  theme(legend.key = element_rect(colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) 
gf
