setwd('/home/piss/Documents/Extreme/R resources/IRM')

library(magrittr) # usefull for pipe %>% operators
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyverse)

library(PissoortThesis) # We load it but we will call the functions by
# the namespace environnement :: to make it more clear.

# Apply the created theme to all ggplots without having to specify it !
theme_set(PissoortThesis::theme_piss())



#################### Uccle ################
###########################################


# From 1833. Free dataset from the KNMI
TX_public <- read.csv("TX Uccle public.csv")
TX_public$TX <- TX_public$TX/10

## From IRM
TXTN_open <- read.csv('Uccle TX TN 1901.csv',sep="")
setnames(TXTN_open,"DAY","Date")

TXTN_closed <- read.csv('Uccle TX TN 1901 closed shelter.csv',sep="")
setnames(TXTN_closed,"DAY","Date")

TXTN_closed$day <- substr(TXTN_closed$Date,7,8)
TXTN_closed$month <- as.numeric(substr(TXTN_closed$Date,5,6))
TXTN_closed$year <- substr(TXTN_closed$Date,1,4)


# Retrieve seasons with our function
# Of course, we are here based on meteorological seasons
TXTN_closed$season <- sapply(TXTN_closed$month,
                             function(x) PissoortThesis::func_season(x))



##### Differences between the public and the IRM datasets ######

# remove missing values and start from 1901 => have same period for both.
TX_public_1901 <- TX_public[TX_public$DATE >= 19010101 & TX_public$Q_TX != 9,]

TXTN_open_compare <- TXTN_open[TXTN_open$Date < 20000101, ]
TXTN_closed_compare <- TXTN_closed[TXTN_closed$Date < 20000101, ]

# Now, we decide to compare the TX
diff_tx_open <- TX_public_1901$TX - TXTN_open_compare$TX
diff_tx_closed <- TX_public_1901$TX - TXTN_closed_compare$TX

library(psych)
describe(diff_tx_closed)
describe(diff_tx_open)

diff_tx_open <- data.frame(difference = diff_tx_open, method = 'open')

diff_tx <- rbind(diff_tx_open,
                 data.frame(difference = diff_tx_closed, method =  "closed"))
ggplot(diff_tx, aes(col = method)) + geom_boxplot(aes(y = difference, x = method)) +
  ggtitle("cccc")

sum(equals(diff_tx_open$diff, 0.0), na.rm = T) / length(diff_tx_open$diff)
sum(equals(diff_tx_closed, 0.0), na.rm = T) / length(diff_tx_closed)
# 53 % of the public data is the same with open shellter, and only 12 % for
# closed shellter. Indeed, this could be problematic (...)
# We go for the closed shellter's data as advised by the IRM.



rm(diff_tx_open, diff_tx_closed, diff_tx_closed, diff_tx_open)

#=======================================================================


# Insert "-" in dat so as they match date values in R
TXTN_closed$Date <- gsub('^(.{4})(.*)$', '\\1-\\2', TXTN_closed$Date)
TXTN_closed$Date <- gsub('^(.{7})(.*)$', '\\1-\\2', TXTN_closed$Date)
TXTN_closed$Date <- as.Date(TXTN_closed$Date)


## NA's in our final dataset ?

sum(is.na(TXTN_closed))   #  0, it's fine !


######################################################################
##################  Analysis in  CLOSED SHELLTER   #################
######################################################################
######################################################################
# From now, all the analysis below are based on temp. recorded with closed shelters.


############## Group temp by month

ggplot(data=TXTN_closed, aes(group = month)) + geom_boxplot(aes(x = month, y = TX)) +
  theme_bw()
# Variance is quite similar accros months
# quite same distribution of the TX accross months

library(RColorBrewer)
# months
ggplot(data = TXTN_closed, aes(TX, colour = as.factor(month))) +
  geom_density(size = 1.1) +  scale_color_discrete()
# !! same smoothing factor for all densities

# seasons
ggplot(data=TXTN_closed,aes(TX,colour=season)) + geom_density(size=1.1)

ggplot(data=TXTN_closed,aes(x=season,y=TX,group=season)) + geom_boxplot()



##########################################################
##################    GEV     ###########################
##########################################################
library(evd)
library(extRemes)


# block length : usual method of 1 year

list_by_years <- split(TXTN_closed, TXTN_closed$year)
# 116 years of data !


##################   MAX  by year.  See ?yearly.extrm()
max_years <- PissoortThesis::yearly.extrm()

####################   MIN
min_years <- PissoortThesis::yearly.extrm(Fun = min, tmp = "TN")



## Get a global view of the serie of extrema
library(xts)
library(dygraphs)
library(zoo)

xtdata0 <- xts(TXTN_closed$TX, order.by = (TXTN_closed$Date), f = 12)
dygraph(xtdata0, main = "(Dynamic) Time series of TX in Uccle",
        xlab = "Date", ylab = "TX") %>% dyRangeSelector()
# Seasonality is remarkable... indications of nonstationarity !...

xtdata <- xts(max_years$df$Max, order.by = as.yearmon(max_years$df$Year), f = 12)
dygraph(xtdata) %>% dyRangeSelector()
# Well another shape when we take


## Plot the yearly maxima together with some "standard" fitting methods

summary(lm1 <- lm(max_years$data ~ max_years$df$Year))
lm1_1 <- lm1$coefficients[1]
lm1_2 <- lm1$coefficients[2]

Broken_lin1 <-  predict(lm(max_years$data[1:75] ~ max_years$df$Year[1:75]) )
Broken_lin2 <-  predict(lm(max_years$data[77:116] ~ max_years$df$Year[77:116]) )

g1 <- ggplot(data = max_years$df,aes(x=Year,y=Max)) + geom_line() +
  geom_point() +
  geom_smooth(method='lm',formula=y~x, aes(colour = "Linear")) +
  geom_line(data = max_years$df[max_years$df$Year %in% 1901:1975,],
            aes(x = Year, colour = "BrokenLinear", y = Broken_lin1),
             size = 1.5, linetype = "twodash") +
  geom_line(data = max_years$df[max_years$df$Year %in% 1977:2016,],
            aes(x = Year, colour = "BrokenLinear", y = Broken_lin2),
            size = 1.5, linetype = "twodash") +
  stat_smooth(method = "loess", se = F, aes(colour = 'LOESS')) +
  labs(title = "Complete Serie of Annual TX in Uccle",
       y = expression( Max~(T~degree*C))) +
  theme(axis.line = element_line(color="#33666C", size = .45)) +
  scale_colour_manual(name="Trend",
                      values=c(Linear="blue", BrokenLinear="cyan", LOESS="red")) +
  theme_piss(legend.position = c(.888, .152)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
g1
# Red line is local polynomial regression fitting curve (loees)
# The (optimal) default method is convenient. See ?loess


ggplot(data = max_years$df,aes(x=Year,y=Max)) + geom_point() +
  geom_abline(intercept = lm1_1,slope=lm1_2,col="red")+ theme_bw()

ggplot(data = max_years$df,aes(x=Year,y=Max)) + geom_point() +theme_bw()
ggplot(data = max_years$df,aes(x=Max)) + geom_histogram() +theme_minimal()
ggplot(data = max_years$df,aes(x=Max)) + geom_density() +theme_bw()



###  what for the minima ??
g2 <- ggplot(data = min_years$df, aes(x=Year,y=Min)) + geom_line() +
  geom_smooth(method='lm',formula=y~x) + geom_point() +
  stat_smooth(method = "loess", se = F, col = 'red' ) +
  labs(title = "Complete Serie of Annual TN in Uccle") +
  theme_piss(20, 15) +
  theme(axis.line = element_line(color="#33666C", size = .45))

# as expected , trend is a bit less strong as for maxima
summary(lm_min <- lm(min_years$data ~ min_years$df$Year))



## Both ?
grid.arrange(g1, g2, nrow = 2)
# and with the full series ? see intro here


###############################
# Please refer to the other code of introduction for the trend analysis with splines
################################



# Preliminaries
###############
gev_tx <- gev.fit(max_years$data)
gev_tx1 <- fgev(max_years$data)
gev_tx2 <- fevd(max_years$data, units = "deg C")
summary(gev_tx2)
distill(gev_tx2)
plot(gev_tx2)

##################
ci(gev_tx2,type="parameter")
ci(gev_tx2,type="parameter",which.par = 3,xrange=c(-.4,.1),method="proflik",verbose=TRUE)
ci(gev_tx2,method="proflik",xrange=c(35,40),verbose=TRUE)
#it is often better to obtain confidence intervals from the profile-likelihood method in
#order to allow for skewed intervals that may better match the sampling df of these parameter
plot(gev_tx2,"trace")

# Gumbel ? ie test if shape parameter is significant
gev_gumb <- fevd(max_years$data,type="Gumbel",units="deg C")
plot(gev_gumb)   # fit is well poorer
plot(gev_gumb,"trace")
lr.test(gev_gumb,gev_tx2)  # p-value=0.001 --> we reject H0 that xi=0

# SHAPE = -0.25 hence Weibull type and the upper endpoint is
upp_endpoint <- gev_tx1$estimate['loc'] -
                gev_tx1$estimate['scale'] / gev_tx1$estimate['shape']
upp_endpoint   ;   max(max_years$data)  # OK greater than the maximum value of the data
# estimated proba of exceeding value >= will be exactly 0. :~ Bound for return level

# simulate data from fitted GEV
sim <- rextRemes(gev_tx2,100)
fit_sim <- fevd(sim)
plot(fit_sim)

######  Manually

# Function to compute the (negative) log-likelihood
'gev.nloglik' <- function(param, dataset){
  mu <- param[1]  ; sigma <- param[2]   ;   xi <- param[3]
  m <- min((1 + (xi * (dataset - mu)/sigma)))   # minimize the negative log-likelihood !!
  if(m < 1e-6) return(as.double(1e6)) # High value for neg.L-Likelihood
  if(sigma < 1e-6) return(as.double(1e6)) # ensure that conditions are fulfilled
  if(xi == 0){    # Gumbel log-likelihood
    loglik <- -length(dataset) * log(sigma) - sum((dataset - mu) / sigma)
     -sum(exp(-((dataset - mu) / sigma)))}
  else{    # Frechet or Weibull log-likelihood
    loglik <- -length(dataset) * log(sigma)
     -(1/xi + 1) * sum(log(1 + (xi * (dataset - mu) / sigma)))
     -sum((1 + (xi * (dataset-mu)/sigma))^(-1/xi))
    }
 return(-loglik)
}

# Or
gev.loglik <- function(theta, data){
  y <- 1 + (theta[1] * (data - theta[2]))/theta[3]
  if((theta[3] < 0) || (min(y) < 0)) {
    ans <- 1e+06
  } else {
    term1 <- length(data) * logb(theta[3])
    term2 <- sum((1 + 1/theta[1]) * logb(y))
    term3 <- sum(y^(-1/theta[1]))
    ans <- term1 + term2 + term3
  }
  ans
}

param <- c(mean(max_years$df$Max),sd(max_years$df$Max),
           0.1 )
# 0.1 starting value for Xi is common (~ mean of all practical applications in meteo), or see Coles

dataset <- max_years$data
nlm(gev.loglik, param, data = dataset,
    hessian=T, iterlim = 1e5)$estimate
#
# fn <- function (param, data)
#   -(gev.nloglik(mu = param[1], sig = param[2], xi = param[3], data ))
#
# param <- c(mean(max_years$df$Max), sd(max_years$df$Max), 0.1 )
# 0.1 starting value for Xi is common (~ mean of all practical applications in meteo)

dataset <- max_years$data
numerical_max <- nlm(gev.nloglik, param, data = dataset)  # Problem ?
numerical_max

nlm(fn, param, data = max_years$data, hessian=T, iterlim = 1e5)
# hessian is exactly the obs. Info Matrix as we dealt with -loglik

optim(param, gev.nloglik, data = dataset, hessian = T, method = "BFGS")

Var_max <- solve(numerical_max$hessian)
sqrt(diag(Var_max))

###############

library(stargazer)

tab <- rbind.data.frame(gev_tx2$results$par, unname(gev_tx1$std.err) )
colnames(tab) <- c("Location", "Scale", "Shape")
rownames(tab) <- c("Estimates", "Std.errors")
tab
stargazer(tab, summary = F)


## Other Methods



#############   Model diagnostics   #################
   ############################################


#################### Probability plot ##################

# get order statistics
max_order <- sort(max_years$data)

#  retrieve the the empirical distribution function
empirical= c()
for(i in 1:length(max_order)){
  empirical[i] <- i/(length(max_years$data)+1)
}
#ecdf(max_data)

# compute Distribution function of the modelled GEV
GEV.DF <- function(data,mu,sigma,xi){
  if(xi == 0)  GEV <- exp(-exp(-((data-mu)/sigma)))
  else  GEV <- exp(-(1+xi*((data-mu)/sigma))^(-1/xi))
 return(GEV)
}

model_est <- c()
for(i in 1:length(max_order)){
  model_est[i] <- GEV.DF(max_order[i],gev_tx1$estimate[1],
                         gev_tx1$estimate[2],gev_tx1$estimate[3])
}
gg_pp <- ggplot(data = data.frame(empirical,model_est),
       aes(x=empirical,y=model_est)) +
  geom_point(shape = 1, col = "#33666C") + geom_abline(intercept=0,slope=1,col="red") +
  theme_piss(16, 11) + labs(y = "Estimated proportions", x = "Normal proportions") +
  ggtitle("PP-plot")

# Fit seems quite well


#################### QQ-plot ##################

# Compute the quantile function (inverse of DF)
model_quantile <- length(max_order)

GEV.INV <- function(data, mu, sigma, xi){
  if(xi==0){  INV <- mu - sigma * log(-log(1 - data))  }
  else{ INV <- mu + (sigma/xi) * (((-log(data))^(-xi))-1)  }
return(INV)
}

for(i in 1:length(max_order)){
 model_quantile[i] <- GEV.INV(empirical[i], gev_tx1$estimate[1],
                              gev_tx1$estimate[2], gev_tx1$estimate[3])
}
gg_qq <- ggplot(data=data.frame(model_quantile,max_order), aes(x=max_order,y=model_quantile)) +
  geom_point(shape = 1, col = "#33666C") + geom_abline(intercept=0,slope=1,col="red") +
  theme_piss(16, 11) + labs(y = "Model quantiles", x = "Normal quantiles") +
  ggtitle("QQ-plot")
# Same conclusion


gridExtra::grid.arrange(gg_qq,gg_pp, nrow = 1)



################## Return Level ######################
#####################################################

# 100-year return level ?
y_100 <- -(1 / log(1 - 1/100))
r_m100 <- gev_tx1$estimate[1] + (gev_tx1$estimate[2] / gev_tx1$estimate[3]) *
  (y_100^gev_tx1$estimate[3] - 1)
# see eq. in the text --> here we take the estimate
GEV.INV(1 - 1/100, gev_tx1$estimate[1], gev_tx1$estimate[2], gev_tx1$estimate[3])
# or directly by our Inverse function actually, at survival probability

return.level(gev_tx2, return.period = c(2, 10, 100, 10000))
## Beware that these assume that data are stationary, while they are not (?) !!
# We remark clearly that these inferences on return levels do not take into account
# the climate warming, ie the fact that mean temperature increases slightly over time
# AND the fact that extremes are more frequent in a climate change [ref] (...)
# Indeed, one shall not expect a TX of 38deg reached only after 10,000 years....

# Standard errors of the estimates
m <- 20  # return period
r_m <- -log(1 - (1/m))
nabla <- matrix(ncol = 1, nrow = 3)
nabla[1,1] <- 1
nabla[2,1] <- -((gev_tx1$estimate[3])^(-1))*(1-(r_m^(-gev_tx1$estimate[3])))
nabla[3,1] <- ((gev_tx1$estimate[2])*
                 ((gev_tx1$estimate[3])^(-2))*(1-((r_m)^(-gev_tx1$estimate[3]))))-
  ((gev_tx1$estimate[2])*((gev_tx1$estimate[3])^(-1))*
     ((r_m)^(-(gev_tx1$estimate[3])))*log(r_m))
sqrt(t(nabla)%*%gev_tx1$var.cov%*%nabla)

# Retrieve "return levels" associated to prob of exceeding thresold from fitted GEV
proba_excess <- pextRemes(gev_tx2, q = c(25, 30, 35, 38, upp_endpoint),
                          lower.tail = F)
(proba_excess)^-1
sum(max_years$data > 35)
# a bit surprizingly, we see that the probability of exce

### Return leel Plot ?
# rlplot(gev_tx)
#
#  rr <- return.level(gev_tx2, return.period = 2:500)
# df_rl <-  data.frame("Return Period" = 2:500,"Return Level" = unname(rr[2:500]))
#
# ggplot(df_rl) + geom_point(aes(x = Return.Period, y = Return.Level),shape = 1)

# See function built in UsedFunc.R
gg_rl <- rl_piss(gev_tx$mle, gev_tx$cov, gev_tx$data)

grid.arrange(pp_gg, gg_rl, nrow = 1)

library(cowplot)

##############################

library(ismev)

gev.diag(gev_tx) # confidence bands seem not too large in the return leel plot !

x <- seq(min(max_years$data), max(max_years$data), length = length(max_years$data))
weib_fit_dens <- evd::dgev(x,loc = gev_tx$mle[1],
                           scale = gev_tx$mle[2], shape = gev_tx$mle[3])

ggplot(data.frame(x,weib_fit_dens)) + geom_density(aes(x = max_years$data)) +
   geom_line(aes(x = x, y = weib_fit_dens), colour = "red") + theme_bw() +
   coord_cartesian(xlim = c(26, 37)) +
   geom_vline(xintercept = min(max_years$data), linetype = 2) +
   geom_vline(xintercept = max(max_years$data), linetype = 2)
# Red is the fitted density while black is the empirical one.

#### Profile likelihood  ( For stationary model !!! )
gev.prof(gev_tx,20,xlow=min(max_data)+7,xup=max(max_data))

gev.profxi(gev_tx,xlow=0,xup=2)

par(mfrow=(c(1,3)))
plot(profile(gev_tx1), ci=c(0.95,0.99))


# Return Levels and empirical Quantiles
gev_tx <- gev.fit(max_years$data)

gg_rl <- PissoortThesis::rl_piss(gev_tx$mle, gev_tx$cov, gev_tx$data)

## Density plots
x <- seq(min(max_years$data)-5, max(max_years$data)+5, length = length(max_years$data))
weib_fit_dens <- evd::dgev(x,loc = gev_tx$mle[1],
                           scale = gev_tx$mle[2], shape = gev_tx$mle[3])

density <- c( "empirical" = "blue", "fitted" = "red")
gg_ds <- ggplot(data.frame(x,weib_fit_dens)) + stat_density(aes(x = max_years$data, col = "empirical"), geom = "line") + ggtitle("Empirical (black) vs fitted (red) density") +
  geom_line(aes(x = x, y = weib_fit_dens, col = "fitted")) + theme_piss(size_p = 5) +
  coord_cartesian(xlim = c(25, 38)) + labs(x = "TX") +
  geom_vline(xintercept = min(max_years$data), linetype = 2) +
  geom_vline(xintercept = max(max_years$data), linetype = 2) + theme_piss(17) +
  scale_colour_manual(name = "Density", values = density) +
  theme(legend.position = c(.888, .82))  +
  guides(colour = guide_legend(override.aes = list(size = 2)))
gg_ds
# Red is the fitted density while black is the empirical one.

gridExtra::grid.arrange(gg_rl, gg_ds, nrow = 1)

##########################################################
# ##################    POT   ########################
##########################################################

library(ismev)

max_all <- TXTN_closed$TX

threshrange.plot(max_all,r=c(28,35))

# mean residual life plot
u <- seq(min(max_all),max(max_all),0.01)
x <- c()
for (i in 1:length(u)) {
  threshold.excess <- max_all[max_all>u[i]]
  x[i] <- mean(threshold.excess-u[i])
}
plot(x~u,type="l",main="MRL plot",ylab="mean excess")

mrlplot(max_all)
mrl.plot(max_all)
mrl.plot(max_all, umin = 28, umax = max(max_all))  # Explore values
# Let's go for 31-32deg as threshold (?), as see ~some linearity in mrl plot

library(evmix)
mrlplot(max_all, tlim = c(20, max(max_all)))
tcplot(max_all,tlim = c(20, max(max_all)) )


### Do not consider the following pack as EVI<0 here.
library(PissoortThesis)
MeanExcess(max_all, plot = T, k =1000)
genHill(max_all)
library(evir)
hill(max_all,start = 10,end = 10000)

# Dispersion Index Plot

diplot(data.frame(obs = seq(1:length(max_all)),max_all), u.range = c(20,36) )
cbind(TXTN_closed$t,TXTN_closed$Date,TXTN_closed$TX)
diplot(cbind(obs=TXTN_closed$t,time=TXTN_closed$Date,TXTN_closed$TX), u.range = c(20,36) )





fit_pot1 <- gpd.fit(max_all, 25)
gpd.diag(fit_pot1)
fit_pot <- fevd(max_all,threshold = 25, type="GP")
plot(fit_pot)
pextRemes(fit_pot,c(25,30,35,36), lower.tail = FALSE)


# Not possible to do POT for whole days in year. max in winter will not be max etc
####################################"

# see distribution of TX wrt "arbitrary' threshold
above_25 <- TXTN_closed[TXTN_closed$TX>25,] #  ( change value of 25)

ggplot(data=as.data.frame(above_25),aes(x=as.factor(month),y=TX)) + geom_boxplot()
# Only months from april to october have TX excessing 25degC


## Clusters ?
above_32 <- TXTN_closed[TXTN_closed$TX>32, ] #  ( change value of 25)
ggplot(above_32) + geom_point(aes(x = Date, y = TX))
# Indeed clusters of exceedance with POT as exceedances are not independent
# ---> stationary ?   One can also remark particular warm years, like 1976
# Declustering techniques ?



####################  Dividing analysis Winter-Summer, seasons, etc...  ###########
###################################################################################


TXTN_closed_s <- TXTN_closed[TXTN_closed$month %in% c(4, 5, 6, 7, 8, 9),]
TXTN_closed_w <- TXTN_closed[TXTN_closed$month %in% c(1, 2, 3, 10, 11, 12),]

describe(TXTN_closed_s)
describe(TXTN_closed_w)

TXTN_closed_s <- data.frame (TXTN_closed_s, period = "1")
TXTN_closed_w <- data.frame (TXTN_closed_w, period = "0")
TXTN_closed_ws <- rbind(TXTN_closed_w, TXTN_closed_s)

ggplot(TXTN_closed_ws, aes(x = TX, col = period)) +
  geom_density() + theme_bw()

gX <- ggplot(TXTN_closed_ws, aes(x = period, y = TX)) + geom_boxplot()
gN <- ggplot(TXTN_closed_ws, aes(x = period, y = TN)) + geom_boxplot()
grid.arrange(gX, gN, nrow = 1)


## Retrieve the pertaining datasets w and s, enabling compute max/min for "seasons" (w or s)
list_years_w <- split(TXTN_closed_w,TXTN_closed_w$year)
list_years_s <- split(TXTN_closed_s,TXTN_closed_s$year)


######    WINTER

max_w <- yearly.extrm(list_years_w)
# Replace the last value which is too low as TX for winter did not occur
max_w$df[116,]$Max <- mean(max_w$df$Max[110:115])

summary(lm_w <- lm(max_w$df$Max ~ max_w$df$Year))

ggplot(data = max_w$df, aes(x = Year,y = Max)) + geom_point() +theme_bw()
ggplot(data = max_w$df, aes(x = Max)) + geom_histogram() +theme_minimal()
ggplot(data = max_w$df, aes(x = Max)) + geom_density() +theme_bw()

#####    SUMMER

max_s <- yearly.extrm(list_years_s)


summary(lm_s <- lm(max_s$df$Max ~ max_s$df$Year))
# As expected, we remark that trend is heavier and "more significant"
# in summer months for TX than in winter months

ggplot(data = max_s$df,aes(x = Year, y = Max)) + geom_point() +theme_bw()
ggplot(data = max_s$df,aes(x = Max)) + geom_histogram() +theme_minimal()
ggplot(data = max_s$df,aes(x = Max)) + geom_density() +theme_bw()



# Comparisons
ggw <- ggplot(data = max_w$df, aes(x = Year, y = Max)) + geom_line(size = 0.3) +
  geom_smooth(method='lm',formula = y~x, size = 0.5) +
  stat_smooth(method = "loess", se = F, col = 'red', size = 0.5) +
  ggtitle("TX For Winter months [10-03]") +
  theme_piss(14, 12)
  theme(legend.title = element_text(colour="#33666C",
                                    size=12, face="bold")) +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.y = element_text(size = 5))

ggs <- ggplot(data = max_s$df ,aes(x = Year,y = Max)) + geom_line(size = 0.3) +
  geom_smooth(method = 'lm',formula = y~x, size = 0.5) +
  stat_smooth(method = "loess", se = F, col = 'red', size = 0.5) +
  ggtitle("TX For Summer months [04-09]") +
  theme_piss(14, 12) +
  theme(legend.title = element_text(colour="#33666C",
                                    size=12, face="bold")) +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.y = element_text(size = 5))
grid.arrange(ggw, ggs, nrow = 2)



###################     For MINIMA     ##############################

## Retrieve the minmimum for datasets w and s
## winter

min_w <- yearly.extrm(list_years_w, Fun = min, tmp = 'TN')


summary(lm_w_min <- lm(min_w$df$Min ~ min_w$df$Year))

ggplot(data =min_w$df ,aes(x = Year, y = Min)) + geom_line() +
  geom_smooth(method = 'lm',formula = y~x) + theme_bw() +
  stat_smooth(method = "loess", se = F, col = 'red')

ggplot(data = min_w$df, aes(x = Year, y = Min)) + geom_point() + theme_bw()
ggplot(data = min_w$df, aes(x = Min)) + geom_histogram() + theme_minimal()
ggplot(data = min_w$df, aes(x = Min)) + geom_density() + theme_bw()

## summer

min_s <- yearly.extrm(list_years_s, Fun = min, tmp = "TN")


summary(lm_s_min <- lm(min_s$df$Min ~ min_s$df$Year))


ggplot(data = min_s$df, aes(x = Year,y = Min)) + geom_point() +theme_bw()
ggplot(data = min_s$df, aes(x = Min)) + geom_histogram() +theme_minimal()
ggplot(data = min_s$df, aes(x = Min)) + geom_density() +theme_bw()


gmin.w <- ggplot(data =min_w$df ,aes(x = Year, y = Min)) + geom_line(size = 0.3) +
  geom_smooth(method = 'lm', formula = y~x, size = 0.5) +
  stat_smooth(method = "loess", se = F, col = 'red', size = 0.5) +
  ggtitle("TN for winter months [10-03]") +
  theme_piss(14, 12) +
  theme(legend.title = element_text(colour="#33666C",
                                    size=12, face="bold")) +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.y = element_text(size = 5))

gmin.s <- ggplot(data = min_s$df, aes(x = Year, y = Min)) + geom_line(size = 0.3) +
  geom_smooth(method = 'lm', formula = y~x, size = 0.5) +
  ggtitle("TN for summer months [04-09]") +
  stat_smooth(method = "loess", se = F, col = 'red', size = 0.5) +
  theme_piss(14, 12) +
  theme(legend.title = element_text(colour="#33666C",
                                    size=9, face="bold")) +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.y = element_text(size = 5))


grid.arrange(gmin.w, gmin.s, nrow = 2)


####All
grid.arrange(ggw, ggs,gmin.w, gmin.s, nrow = 4)


var(max_s$data)  ;  var(max_w$data)
var (min_s$data)  ;  var(min_w$data)
# Variance well more for minimi in winter than in summer months.


# Or more convenient to work with NEGATED data for TN ?
neg_min_s <- min_s$df   ;   neg_min_s$negMin <- -neg_min_s$Min
neg_min_w <- min_w$df   ;   neg_min_w$negMin <- -neg_min_w$Min

gnegmin.w <- ggplot(data = neg_min_w ,aes(x = Year, y = negMin)) + geom_line() +
  geom_smooth(method = 'lm', formula = y~x) + theme_bw() +
  stat_smooth(method = "loess", se = F, col = 'red') + ggtitle("Min for winter months")
gnegmin.s <- ggplot(data = neg_min_s, aes(x = Year, y = negMin)) + geom_line() +
  geom_smooth(method = 'lm', formula = y~x) + theme_bw() + ggtitle("Min for summer months") +
  stat_smooth(method = "loess", se = F, col = 'red')
grid.arrange(gnegmin.w, gnegmin.s, nrow = 2)


#################
# As expected, the trend is (slightly) heavier for TX in warm months
#  and heavier for TN in cold month. Test if it is really significant

# http://stats.stackexchange.com/questions/33013/what-test-can-i-use-to-compare-slopes-from-two-or-more-regression-models?rq=1
summary(glm(TXTN_closed_ws$TX ~ max_s$df$Year * TXTN_closed_ws$period))

summary(lm(TXTN_closed_ws$TX~TXTN_closed_ws$period * TXTN_closed_ws$year))

summary(lm(TXTN_closed_ws$TX~TXTN_closed_ws$period))

#



###############################################################################
##################    POT taken on divided datasets   ########################
#############################################################################
library(fExtremes)

TX_all_s <- TXTN_closed_s$TX
TX_all_w <- TXTN_closed_w$TX

############# max TX on summer dataset

threshrange.plot(TX_all_s, r = c(25,35))

# mean residual life plot
u <- seq(min(TX_all_s),max(TX_all_s),0.01)
x <- c()
for (i in 1:length(u)) {
  threshold.excess <- TX_all_s[TX_all_s > u[i] ]
  x[i] <- mean(threshold.excess - u[i])
}
plot(x~u,type="l",main="MRL plot",ylab="mean excess")
mrl.plot(TX_all_s)  # Decreasing : we have LTE distribution (xi negative)
# Let's go for 32-34deg as threshold, as see ~some linearity in mrl plot

above_thres_s <- TX_all_s[TX_all_s>32]

fit_pot_s_1 <- gpd.fit(TX_all_s,32)
gpd.diag(fit_pot_s_1)
fit_pot_s_11 <- fevd(TX_all_s,threshold = 32,type="GP")
plot(fit_pot_s_11)
pextRemes(fit_pot,c(25,30,35,36),lower.tail = FALSE)



par(mfrow=c(2,1))
acf(above_thres_s,lag.max = length(above_thres_s))
pacf(above_thres_s,lag.max = length(above_thres_s))
# This clearly increase when we decrease the threshold...

plot(above_thres_s[1:length(above_thres_s)-1],above_thres_s[2:length(above_thres_s)])
# Very slight sign of dependence....




save.image("data.Rdata") #  To load in the other scripts



