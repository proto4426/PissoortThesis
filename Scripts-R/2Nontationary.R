#setwd('/home/piss/Documents/Extreme/R resources/IRM')
#setwd("C:\\Users\\Piss\\Documents\\LINUX\\Documents\\Extreme\\R resources\\IRM")
load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1.Rdata")
setwd('/home/proto4426/Documents/Extreme/R resources/IRM')
library(tidyverse)
library(gridExtra)

library(evd)
library(extRemes)
library(ismev)
#library(extremeStat)

library(PissoortThesis)


################### Introduction :
##################################

######### Division and analysis of the 4 real seasons separately ? #######
                                      ##########################
#### For TX

ggplot(TXTN_closed[1:2000, ]) + geom_point(aes(x = Date, y = TX))  # change value
# As we have seen, strong reccurent pattern occurs due to seasons.



TXTN_cl_wint <- TXTN_closed[TXTN_closed$season == "Winter", ]
TXTN_cl_spring <- TXTN_closed[TXTN_closed$season == "Spring", ]
TXTN_cl_summ <- TXTN_closed[TXTN_closed$season == "Summer", ]
TXTN_cl_autu <- TXTN_closed[TXTN_closed$season == "Autumn", ]

# Add column which allows to retrieve the index for plotting
TXTN_cl_wint <- cbind(TXTN_cl_wint, index = 1:nrow(TXTN_cl_wint))
TXTN_cl_spring <- cbind(TXTN_cl_spring, index = 1:nrow(TXTN_cl_spring))
TXTN_cl_summ <- cbind(TXTN_cl_summ, index = 1:nrow(TXTN_cl_summ))
TXTN_cl_autu <- cbind(TXTN_cl_autu, index = 1:nrow(TXTN_cl_autu))

ga <- ggplot(TXTN_cl_autu[1:1000,]) + geom_point(aes(x = index, y = TX)) +
  ggtitle("TX For autmun")
gw <- ggplot(TXTN_cl_wint[1:1000,]) + geom_point(aes(x = index, y = TX)) +
  ggtitle("TX for winter")
gs <- ggplot(TXTN_cl_spring[1:1000,]) + geom_point(aes(x = index, y = TX)) +
  ggtitle("TX for spring")
gsu <- ggplot(TXTN_cl_summ[1:1000,]) + geom_point(aes(x = index, y = TX)) +
  ggtitle("TX for summer")
# (change)
grid.arrange(ga, gw, gs, gsu)
# Not sufficient to model seasons separately


# Or for the whole series
library(xts)
library(dygraphs)
library(zoo)

xtdataw <- xts(TXTN_cl_wint$TX, order.by = as.Date(TXTN_cl_wint$index), f = 12)
dygraph(xtdataw) %>% dyRangeSelector()
xtdatas <- xts(TXTN_cl_spring$TX, order.by = as.Date(TXTN_cl_spring$index), f = 12)
dygraph(xtdatas) %>% dyRangeSelector()
xtdatasu <- xts(TXTN_cl_summ$TX, order.by = as.Date(TXTN_cl_summ$index), f = 12)
dygraph(xtdatasu) %>% dyRangeSelector()
xtdataau <- xts(TXTN_cl_autu$TX, order.by = as.Date(TXTN_cl_autu$index), f = 12)
dygraph(xtdataau) %>% dyRangeSelector()
# We easily remark that, even when divided by season, signs of
# nonstationarities remain. This is of course well less strong than before.




######### Warmest months for TX (july and augustus) #############
         ######################

txtn_cl_warm <- TXTN_cl_summ[TXTN_cl_summ$month %in% c(7,8), -8 ]
txtn_cl_warm <- cbind(txtn_cl_warm, index = 1:nrow(txtn_cl_warm))

gw <- ggplot(txtn_cl_warm[1:2000,]) + geom_line(aes( x = index, y = TX)) +
  ggtitle("TX for july and august") +
  geom_hline(yintercept = mean(txtn_cl_warm$TX), col= "red")
gw


list_years_w <- split(txtn_cl_warm, txtn_cl_warm$year)
yearl_ext_warm <- yearly.extrm(list_years_w)
summary( gev_warm <- fevd(yearl_ext_warm$data, units = "deg C") )



######### Coolest months for TN (Januari and februari) #############
         ######################

txtn_cl_cool <- TXTN_cl_wint[TXTN_cl_wint$month %in% c(1,2), -8 ]
txtn_cl_cool <- cbind(txtn_cl_cool, index = 1:nrow(txtn_cl_cool))

gc <- ggplot(txtn_cl_cool[1:2000, ]) + geom_line(aes( x = index, y = TN)) +
  ggtitle("TN for januari and februari") +
  geom_hline(yintercept = mean(txtn_cl_cool$TN), col= "red")
gc

### Change sign of the values and work with maxima !
tn_cl_cool_neg <- txtn_cl_cool
tn_cl_cool_neg$TN <- -(tn_cl_cool_neg$TN)
list_years_c <- split(tn_cl_cool_neg, tn_cl_cool_neg$year)
yearl_ext_cool <- yearly.extrm(list_years_c, Fun = min, tmp = "TN")
summary( gev_cool <- fevd(yearl_ext_cool$data, units = "deg C") )


#####
grid.arrange(gw, gc, nrow = 2)


# Here, we remark with this method that assumption of stationarity is much more
# fulfilled. However, we must also think about the loss of information ....
mrl.plot(txtn_cl_warm$TX)



##############################################################################
#####################    Stationary Series    ################################
############   Relaxing the independence assumption   ########################
##############################################################################


### GEV : Detect autocorrelation/dependence

PissoortThesis::acfpacf_plot(max_years$data, mfrow = c(1, 2),
                            main1 = "ACF for the Annual TX (block maxima)",
                            main2 = "PACF for the Annual TX (block maxima)")
# We can detect/suspect the seasonality with the ~harmonic oscillation
# For block maxima,as expected, autocorrelation is weak (but present!)

plot(max_years$data[1:length(max_years$data) - 1],
     max_years$data[2:length(max_years$data)])
# Seems really independent. For GEV, it seems okay.
# See it for summer months only :
PissoortThesis::acfpacf_plot(max_s$data)


##### See POT in other code ---> clusters
PissoortThesis::acfpacf_plot(above_25$TX)

extremalindex(max_all, threshold = 32) # Method interval was from Segers !!!!!!
clust <- extremalindex(max_all, threshold = 30)
# obviously, dependence increase (theta decreases) as threshold as more extremes


# Controlling number of cluster max,
fpot(max_all, model= "pp",threshold = 32, cmax = T, r = 0,
      start = list(loc = 30, scale = 1, shape = 0), npp = 365.25 )


#( see extReme 2.0 pdf )
atdf(max_all,0.99) # O.95 is u in F^-1(u) over which we compute atdf
# u close to 1 but high enough to incorporate enough data.


## Clusters ?
above_30 <- TXTN_closed[TXTN_closed$TX>30, ] #  ( change value of 30)

ggplot(above_30, aes(x = Date, y = TX)) +
  geom_point() + theme_bw() +
  labs(title = "Occurences of TX above 30°c",
       x = "Date (year)", y= "TX (°c)") +
  theme(plot.title = element_text(size = 22, hjust=0.5,
                                  colour = "#33666C", face="bold")) +
  theme(axis.title = element_text(face = "bold", size= 17,
                                  colour = "#33666C")) +
  theme(axis.line = element_line(color="#33666C", size = .45))  +
  geom_vline(aes(xintercept = as.numeric(as.Date("1911-07-27"))),
             col = "red", size = 0.15) +
  geom_vline(aes(xintercept = as.numeric(as.Date("1976-06-30"))),
             col = "red", size = 0.15)



############################################################
############## Dealing with Nonstationnarity ###############
############################################################

TXTN_closed$Date

lmm <- lm(TXTN_closed$TX ~ TXTN_closed$Date * TXTN_closed$season)
lm0 <- lm(TXTN_closed$TX ~ TXTN_closed$season )
lm1 <- lm(TXTN_closed$TX ~ TXTN_closed$season + TXTN_closed$Date)
summary(lmm)   ;   summary(lm0)  ;     summary(lm1)

df <- data.frame(obs = TXTN_closed$TX[1:1500], season = predict(lm0)[1:1500],
                 season_trend = predict(lm1)[1:1500], interac = predict(lmm)[1:1500],
                 Date = TXTN_closed$Date[1:1500])

ggplot(df) + geom_point(aes(x = Date, y = obs)) +
  geom_line(aes(x = Date, y = season, colour = "red")) +
  geom_line(aes( x = Date, y = season_trend, colour = "blue")) +
  geom_line(aes(x = Date, y = interac, colour = "green"))


###############   GEV   ###########################

# We have seen this reccurent pattern... First, simple linear trend ?
# As we have seen with the LR on all the TX, we expect this to be
#significatie here too (?!)

ti <- matrix (ncol = 1, nrow = length(max_years$data))
#ti[ ,1] <- rep(1, length(max_years$data))
ti[ ,1] <- seq(1, length(max_years$data),1)
gev_nonstatio <- ismev::gev.fit(max_years$data, ydat = ti , mul = c(1))
# Value of b_1 is ~~the same than this obtained with linea regression

# Comparing it with likelihood of stationary model
gev_statio <- gev.fit(max_years$data)
print(khi.obs <- 2 *( (-gev_nonstatio$nllh[1]) - (-gev_statio$nllh[1]) ))
qchisq(.05, df = 1, lower.tail = F)
# It is significant.  P-val ?
pchisq(khi.obs, df = 1, lower.tail = F)
# Not exactly the same as with LR, but same result

## More complex model ?  Quadratic model
ti2 <- matrix(ncol = 2, nrow = length(max_years$data))
#ti2[ ,1] <- rep(1, length(max_years$data))
ti2[ ,1] <- seq(1, length(max_years$data), 1)
ti2[ ,2] <- (ti2[,1])^2
gev_nonstatio2 <- ismev::gev.fit(max_years$data, ydat = ti2, mul = c(1, 2))
pchisq(2 *( (-gev_nonstatio2$nllh[1]) - (-gev_nonstatio$nllh[1]) ), df = 1,
       lower.tail = F)
# compared with the khi-2 as above, it is not statisticaly necessary

## 'Cubic' model ?
ti3 <- matrix(ncol = 3, nrow = length(max_years$data))
ti3[ ,1] <- seq(1, length(max_years$data), 1)
ti3[ ,2] <- (ti3[,1])^2
ti3[ ,3] <- (ti3[,1])^3
gev_nonstatio3 <- gev.fit(max_years$data/1000, ydat = ti3, mul = c(1, 2, 3))
gev_nonstatio3 <- PissoortThesis::gev.fit2(max_years$data, ydat = ti3,
                                           mul = c(1, 2, 3),
                                           browser = T, solve.tol = 1e-25)
# System is singular if we do not scale the data.  But if we scale, the are still
#NaNs produced in the covariance matrix. We needed to change the tolerance of the
#solve() function for the underlying hessian (covariance) matrix.

pchisq(2 *( (-gev_nonstatio3$nllh[1]) - (-gev_nonstatio$nllh[1]) ),
        df = 2, lower.tail = F)


## Linear trend with varying scale parameter
ti_sig <- ti
gev_nonstatio_sig <- ismev::gev.fit(max_years$data, ydat = ti_sig,
                                    mul = c(1), sigl = c(1))
pchisq(2 *( (-gev_nonstatio_sig$nllh[1]) - (-gev_nonstatio$nllh[1]) ),
       df = 1, lower.tail = F)
# No reason to vary scale with time

### Nonstationary scale parameter ?
ti_sig <- matrix(ncol = 1, nrow = length(max_years$data))
ti_sig[,1] <- seq(1, length(max_years$data), 1)
gev_nstatio_scale <- ismev::gev.fit(max_years$data, ydat = ti_sig,
                                    sigl = 1, siglink = exp)
pchisq(2 *( (-gev_nstatio_scale$nllh[1]) - (-gev_statio$nllh[1]) ), df = 1,
       lower.tail = F)
# This is not significant



##### Diagnostics
# ismev by default allows for nonstationary model in gev.diag()
gev.diag(gev_nonstatio)
# problem for high quantiles in the qq-plot (as "usual")

# Produce the diagnostic plots in ggplot2
Empirical <- c()
for(i in 1:length(max_order)){
  Empirical[i] <- i/(length(max_years$data)+1)
}
gg_pp.trans <- ggplot(data = data.frame(Empirical,
                                  Model = exp(-exp(-sort(gev_nonstatio$data)))),
                aes(x = Empirical, y = Model)) +
  geom_point(shape = 1, col = "#33666C") +
  geom_abline(intercept=0,slope=1,col="red") +
  theme_piss(16, 11) +
  labs(y = "Estimated proportions", x = "Empirical proportions") +
  ggtitle("Residual PP-plot")
gg_qq.trans <- ggplot(data = data.frame(Model = sort(gev_nonstatio$data),
                                  Empirical = -log(-log(Empirical))),
                aes(x = Empirical, y = Model)) +
  geom_point(shape = 1, col = "#33666C") +
  geom_abline(intercept=0,slope=1,col="red") +
  theme_piss(16, 11) +
  labs(x = "Model quantile", y = "Empirical quantile") +
  ggtitle("Residual QQ-Plot (Gumbel Scale)")
gridExtra::grid.arrange(gg_pp.trans, gg_qq.trans, nrow = 1)




#####   Return levels for model with linear trend. m = 10
       ##############

rl_10_lin <- PissoortThesis::return.lvl.nstatio(max_years$df$Year,
                                               gev_nonstatio, t = 257, m = 25)
gg_rlAll <- rl_10_lin$g +
  geom_vline(xintercept = 2066, linetype = 2, col = 1, size = 0.15) +
  geom_vline(xintercept = 2216, linetype = 2, col = 1, size = 0.15) +
  geom_vline(xintercept = 2274, linetype = 2, col = 1, size = 0.15) +
  geom_vline(xintercept = ( max(max_years$df$Year) + length(max_years$data) ),
             linetype = 2, col = 2) +
  labs(title = "Return levels with return period of 25 years",
       y = "25-years Return Level",
       x = "Year \n (prediction horizon)") +
  theme_piss(size_p = 12, size_c = 10) +
  scale_x_continuous(breaks = c(2016, 2066, 2132, 2216, 2016+length(rl_10_lin$rl) ),
                labels = c("2016 \n (0)","2066 \n (50)", "2132 \n (116)" , "2216 \n (200)",
                           paste(as.character(2016+length(rl_10_lin$rl)), "\n (257)") ))#+ coord_cartesian(ylim = c())
gg_rlAll
# Note the Increase of the return level with time (due to trend)
# Avec le trend qu'on a pour l'instant, dans environ 300 ans on depassera
# la valeur de 50degres maximum tout les 10ans en moyenne a uccle....
rl_10_lin$rl
#caution should be exercised in practice concerning whether or not it is believable
#for the upward linear trend in maximum temperatures to continue to be valid.
rl_10_lin$rl[116] ;  max(TXTN_closed$TX) # Very close to the maximum of the series


# And if we start at years ~ 65 (donnes plus fiables, + debut du RC)
ti_65 <- matrix (ncol = 1, nrow = 51)
ti_65[ ,1] <- seq(66, 116, 1)
gev_nstatio_65 <- gev.fit(max_years$data[66:116], ydat = ti_65 , mul = 1)
rl_10_lin65 <- PissoortThesis::return.lvl.nstatio(max_years$df$Year[66:116],
                                                 gev_nstatio_65, t = 250, m = 25)
gg_rlSmall <- rl_10_lin65$g +
  geom_vline(xintercept = 2016 + length(max_years$data), linetype = 2, col = 2) +
  labs(title = "Return levels with return period of 25 years") +
  theme_piss(size_p = 13)
# We see that the trend is much heavier, and hence return lvls are too...

grid.arrange(gg_rlAll, gg_rlSmall, nrow = 1)




###########  POT ########################
########################################


### Seasonal model ?

t <- 1:length(TXTN_closed$TX)
seas_effect <- sin(2 * pi * t / 365.25)

summary(lm_seas_sin <- lm(TXTN_closed$TX ~ seas_effect))

df_seas_sin <- cbind(TXTN_closed, seas_effect, t, pred = predict(lm_seas_sin))
ggplot(df_seas_sin[1:1500,]) + geom_point(aes(x = t, y = TX)) +
  geom_line(aes(x = t, y = seas_effect)) + geom_line(aes(x = t, y = pred))


nls.mod <- nls(max_all  ~ a + b * cos(2*pi*t/365.25 - c),
               start = list(a = 1, b = 1, c = 1))
co <- coef(nls.mod)
f <- function(x, a, b, c) {a + b*cos(2*pi*x/365.25 - c) }
ggplot(cbind.data.frame(TXTN_closed[1:5000, ], pred = predict(nls.mod)[1:5000])) +
  geom_point(aes(x = Date, y = TX), size = 0.4) + geom_line(aes(x = Date, y = pred), col = "red")



#### Varying the threshold !

u <- predict(nls.mod)  # (see seasonal model above)
gpd_varu <- fevd(TX, TXTN_closed, threshold = u, verbose = T,
                 type = "GP", time.units = "365/year", units = "deg C")

above_u <- TXTN_closed
above_u$exceed <- above_u$TX - u
above_u1 <- above_u[above_u$exceed>0,]
ggplot(above_u1[1:3000, ], aes(x = Date, y = TX)) + geom_point()
# Lots of exceedances (>2000) but they seem independent

# make the threshold heavier. ( Study bias-variance tradeoff wrt u for the choice)
nls.mod1 <- nls(max_all + 10 ~ a + b * cos(2*pi*t/365.25 - c),
               start = list(a = 1, b = 1, c = 1))
above_u$exceed_red <-  above_u$TX - predict(nls.mod1)
above_ured <- above_u[above_u$exceed_red > 0, ]

ggplot(cbind.data.frame(TXTN_closed[1:5000, ], pred = predict(nls.mod1)[1:5000])) +
  geom_point(aes(x = Date, y = TX), size = 0.4) + geom_line(aes(x = Date, y = pred), col = "red")
ggplot(above_ured, aes(x = Date, y = TX)) + geom_point() +
  geom_smooth(method = "lm",formula = y~x)
table(above_ured$year)   # More values in the end, but more strong excess in the start
# must decluster this

gpd_varu_red <- fevd(TX, TXTN_closed, threshold = predict(nls.mod1),
                     verbose = T, type = "GP", time.units = "365/year", units = "deg C")
plot(gpd_varu_red)
gpd_varu_red <- gpd.fit(max_all, predict(nls.mod1))
gpd.diag(gpd_varu_red)


#######  Declustering

extremalindex(max_all, threshold = predict(nls.mod1), method = "intervals")
# Far from the independent serie....

# Here, value for threshold must be fixed....
fpot(max_all, model= "pp", threshold = predict(nls.mod1), cmax = T, r = 0,
     start = list(loc = 30, scale = 1, shape = 0), npp = 365.25 )

unique(decl_varu <- clusters(max_all,u = predict(nls.mod1), cmax = T, rlow = 0))

ggplot(data.frame(index = names(decl_varu), TX = unname(decl_varu)),
       aes(x = index, y = TX)) + geom_point()
# The series seems much more interesting now !




######

gpd_linearmu <- gpd.fit(max_all,threshold = 32, ydat = ti , mul = 1)

gpd_simple <- fevd(max_all, threshold = 32, type = "GP",
                   time.units = "365/year", units = "deg C")
gpd_seas <- fevd(max_all, threshold = 32,    # change
                type = "GP", time.units = "365/year", units = "deg C")
lr.test(gpd_simple, gpd_seas)

## Trend ?
summary(lm(above_30$TX ~ above_30$t))
ggplot(above_30,aes(x = Date, y = TX)) + geom_point() +
  geom_smooth(method='lm',formula=y~x)
# not really for the exceedances ....



##################  Point Process  ############################
###############################################################


gev_pp <- fevd(max_years$data, units = "deg C", threshold = 32,
               type = "PP", verbose = T)
plot(gev_pp)
# Beware : For non-stationary models, the density plot is not provided.
plot(gev_pp, "trace")
ci(gev_pp, type = "parameter")


threshrange.plot(max_all, r = c(25, 34), nint = 20, type = "PP")
# 32 semms a good choice

# trend ?
attach(TXTN_closed)
gev_pp <- fevd(Max, max_years$df,  threshold = 32,
               type = "PP", verbose = T, threshold.fun = I(Year),
               location.fun = I(Year), time.units = "365/year")
detach(TXTN_closed)


TXTN_closed$t <- t

gev_pp1 <- fevd(TX, TXTN_closed, threshold = 32, location.fun =
                ~ cos(2 * pi * t / 365.25) + sin(2 * pi * t / 365.25),
               type = "PP", units = "deg C")
lr.test(gev_pp, gev_pp1)



pp_coles <- fpot(max_all, model= "pp",threshold = 32,
                 start = list(loc = 30, scale = 1, shape = 0),
                 npp = 365.25 )
plot(pp_coles)
excess_coles <- fpot (max_all, threshold = 32, npp = 365.25)
# Same results for the shape, but no more location parameter
# since it conditions on the exceedances. Different scale



#save.image("data1.Rdata")
