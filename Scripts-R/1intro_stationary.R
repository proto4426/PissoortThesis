#setwd('/home/piss/Documents/Extreme/R resources/IRM')
#setwd("C:\\Users\\Piss\\Documents\\LINUX\\Documents\\Extreme\\R resources\\IRM")
#load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1.Rdata")

setwd('/home/proto4426/Documents/Extreme/R resources/IRM')

library(data.table)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

library(PissoortThesis) # We load it but we will call the functions by
# the namespace environnement :: to make it more clear.

# Apply the created theme to all ggplots without having to specify it !
theme_set(PissoortThesis::theme_piss())



#################### Read Data ###########
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



######################################################################
##################  Analysis in  CLOSED SHELLTER   #################
######################################################################
# From now, all the analysis below are based on temp. recorded with closed shelters.


# Insert "-" in dat so as they match date values in R
TXTN_closed$Date <- gsub('^(.{4})(.*)$', '\\1-\\2', TXTN_closed$Date)
TXTN_closed$Date <- gsub('^(.{7})(.*)$', '\\1-\\2', TXTN_closed$Date)
TXTN_closed$Date <- as.Date(TXTN_closed$Date)


# NA's in our final dataset ?
sum(is.na(TXTN_closed))   #  0, it's fine !


############## Group temperatures by month and inspect the distributions

ggplot(data=TXTN_closed, aes(group = month)) +
  geom_boxplot(aes(x = month, y = TX)) +
  theme_bw()
# Variance is quite similar accros months
# quite same distribution of the TX accross months

## Violin-plots
dodge <- position_dodge(width = 0.4)
gv1 <- ggplot(TXTN_closed,aes(x = season, y = TX)) +
  geom_jitter(color='red', size = .6, alpha=0.99,width = 0.2) +
  geom_violin(fill = "lightseagreen", alpha=0.7, draw_quantiles = T,
              position = dodge, width = 1.8) +
  geom_boxplot(width=.06, outlier.colour=NA, position = dodge) +
  labs(title = 'Violin-plots for daily Max. t°c by seasons',
       x = "Season", y = expression( Maximum~T~degree*C)) +
  theme_piss( size_p = 16)

## Density plots
ggplot(data = TXTN_closed, aes(TX, colour = as.factor(month))) +
  geom_density(size = 1.1) +
  scale_color_discrete() +
  theme_bw()
# !! same smoothing factor for all densities

summer <- TXTN_closed[TXTN_closed$season == "Summer", ]
spring <- TXTN_closed[TXTN_closed$season == "Spring", ]
winter <- TXTN_closed[TXTN_closed$season == "Winter", ]
autumn <- TXTN_closed[TXTN_closed$season == "Autumn", ]

m_summer <- mean(summer$TX)
m_spring =  mean(spring$TX)
m_winter <- mean(winter$TX)
m_autum <- mean(autumn$TX)

gd1 <- ggplot(data=TXTN_closed, aes(TX,fill = season, colour=season)) +
  geom_density(alpha = .1, size=1.1) +
  scale_fill_brewer(palette = "Set1" )+
  scale_color_brewer(palette= "Set1") +
  geom_hline(yintercept=0, colour="white", size=1.1) +
  labs(title = 'Kernel Densities for daily Max. t°c by seasons',
       y = "Density", x = expression( Maximum~T~degree*C)) +
  theme_piss(legend.position = c(0.9, .82), size_p = 16) +
  geom_vline(xintercept = m_summer, colour = "darkgreen", linetype = 2) +
  geom_vline(xintercept = m_spring, colour = "blue", linetype = 2) +
  geom_vline(xintercept = m_winter, colour = "violet", linetype = 2) +
  geom_vline(xintercept = m_autum, colour = "red", linetype = 2)


## violin and density plots together
grid.arrange(gv1, gd1, nrow = 1)





# block length : the usual method is one 1 year
list_by_years <- split(TXTN_closed, TXTN_closed$year)
# Then, we have 116 years of data ! (1901 to 2016)


## Retrieve the max (TN) in each year
max_years <- PissoortThesis::yearly.extrm()

##  min (TN)
min_years <- PissoortThesis::yearly.extrm(Fun = min, tmp = "TN")


######  Get a global and dynamic view of the serie of extremes
library(xts)
library(dygraphs)
library(zoo)

xtdata0 <- xts(TXTN_closed$TX, order.by = (TXTN_closed$Date), f = 12)
dygraph(xtdata0, main = "(Dynamic) Time series of TX in Uccle",
        xlab = "Date", ylab = "TX") %>% dyRangeSelector()
# Seasonality is remarkable... indications of nonstationarity !...  (to see later)
xtdata <- xts(max_years$df$Max, order.by = as.yearmon(max_years$df$Year), f = 12)
dygraph(xtdata) %>% dyRangeSelector()


## Plot the yearly maxima together with some "usual" fitting methods :
#linear regression, LOESS and broken linear regression
lm1 <- lm(max_years$data ~ max_years$df$Year)
lm1_1 <- lm1$coefficients[1]
lm1_2 <- lm1$coefficients[2]
Broken_lin1 <-  predict(lm(max_years$data[1:75] ~ max_years$df$Year[1:75]) )
Broken_lin2 <-  predict(lm(max_years$data[77:116] ~ max_years$df$Year[77:116]) )

gg_trends <- ggplot(data = max_years$df,aes(x=Year,y=Max)) +
  geom_line() + geom_point() +
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
  scale_colour_manual(name="Trend", values=c(Linear="blue",
                                             BrokenLinear="cyan",
                                             LOESS="red")) +
  theme_piss(legend.position = c(.92, .12)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
gg_trends
# Red line is local polynomial regression fitting curve (loees)
# The (optimal) default method is convenient. See ?loess

# Check shaded grey area around linear trend is indeed 95% pointwise ci (in prediction)
summary(lm1)
# x <- seq(1901, 2016, by = 1e-2)
# predict(lm1, x[sample(x, 116)])
# predict(lm1, max_years$df$Year)

conf_int <- predict(lm1, data.frame(max_years$df$Year),
                    level = 0.95, interval="predict", se.fit = T,
                    simultaneous = F)
conf_int$fit
EnvStats::pointwise(conf_int, coverage = 0.95, simultaneous = TRUE)

# These are not the same intervals ! The interval from the shaded grey area are actually
# + - 1.96 times standard errors. We check it :
df_conf <- data.frame(conf_int$fit, Year = max_years$df$Year)
gg_trends_conf <- gg_trends +
  geom_line(aes(x = Year, y = lwr), data = df_conf,
            linetype = 2 , col = "blue") +
  geom_line(aes(x = Year, y = upr), data = df_conf,
            linetype = 2 , col = "blue")
df_conf_se <- data.frame(se = conf_int$se.fit, df_conf)
gg_trends_conf +
  geom_line(aes(x = Year, y = fit + 1.96*se), data = df_conf_se,
            linetype = 2 , col = "orange") +
  geom_line(aes(x = Year, y = fit - 1.96*se), data = df_conf_se,
            linetype = 2 , col = "orange")
# Indeed, our thoughts are verified.
# pointwise ci at year 1920 :
df_conf_se$fit[1920-1900] + c(-1.96 * df_conf_se$se[1920-1900],
                              1.96 * df_conf_se$se[1920-1900] )
#  year 2016 :
df_conf_se$fit[2016-1900] + c(-1.96 * df_conf_se$se[2016-1900],
                              1.96 * df_conf_se$se[2016-1900] )



# Histogram and density plots of the annual maxima
ggplot(data = max_years$df, aes(x=Max)) +
  geom_histogram() +
  theme_minimal()
ggplot(data = max_years$df, aes(x=Max)) +
  geom_density() +
  theme_bw()




#### What for the minima ??
g2 <- ggplot(data = min_years$df, aes(x=Year,y=Min)) +
  geom_line() +
  geom_smooth(method='lm',formula=y~x) +
  geom_point(col = "blue") +
  stat_smooth(method = "loess", se = F, col = 'red' ) +
  labs(title = "Complete Serie of Annual TN in Uccle") +
  theme_piss(20, 15) +
  theme(axis.line = element_line(color="#33666C", size = .45))

# as expected , trend is a bit less strong as for maxima :
summary(lm_min <- lm(min_years$data ~ min_years$df$Year))


## Both maxima and minima on one plot to compare :
grid.arrange(gg_trends, g2, nrow = 2)




###############################
# Please refer now to the other code of introduction for the trend analysis
# with GAM and splines, called intro_trends(splines).R
################################



##########################################################
##################  Analysis by  GEV  ####################
##########################################################
library(evd)
library(extRemes)
library(ismev)


# Fitting GEV by MLE
####################

gev_tx <- gev.fit(max_years$data)
gev_tx1 <- fgev(max_years$data)
gev_tx2 <- fevd(max_years$data, units = "deg C")
# Same results from the different packages : OK
summary(gev_tx2)
distill(gev_tx2)
plot(gev_tx2)

##################
ci(gev_tx2,type="parameter")
ci(gev_tx2,type="parameter",which.par = 3,xrange=c(-.4,.1),method="proflik",verbose=TRUE)
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
sim <- rextRemes(gev_tx2, 100)
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
# Or, more straightforward for the positive log-likelihood
'gev.loglik' <- function(theta, data){
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
# 0.1 starting value for Xi is usual (~ mean of all practical applications in meteo),
#or see Coles (2001)

dataset <- max_years$data
nlm(gev.loglik, param, data = dataset,
    hessian=T, iterlim = 1e5)$estimate
# fn <- function (param, data)
#   -(gev.nloglik(mu = param[1], sig = param[2], xi = param[3], data ))
#
# param <- c(mean(max_years$df$Max), sd(max_years$df$Max), 0.1 )


dataset <- max_years$data
numerical_max <- nlm(gev.nloglik, param, data = dataset)  # Problem ?
numerical_max

nlm(fn, param, data = max_years$data, hessian=T, iterlim = 1e5)
# hessian is exactly the obs. Info Matrix as we dealt with -loglik

optim(param, gev.nloglik, data = dataset, hessian = T, method = "BFGS")

Var_max <- solve(numerical_max$hessian)
sqrt(diag(Var_max))


## Profile likelihood
gev.profxi(gev_tx,xlow=0,xup=2)

par(mfrow=(c(1,3)))
plot(profile(gev_tx1), ci=c(0.95,0.99))

###############
library(stargazer)
tab <- rbind.data.frame(gev_tx2$results$par, unname(gev_tx1$std.err) )
colnames(tab) <- c("Location", "Scale", "Shape")
rownames(tab) <- c("Estimates", "Std.errors")
stargazer(tab, summary = F)


### Other Method of estimation : PWM estimator

library(fExtremes)
fit_pwm <- gevFit(max_years$data,type = "pwm" )

library(lmomco)  # Quite big package
# pwm.gev(max_years$data)


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
  geom_point(shape = 1, col = "#33666C") +
  geom_abline(intercept=0,slope=1,col="red") +
  theme_piss(16, 11) +
  labs(y = "Estimated proportions", x = "Empirical proportions") +
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
gg_qq <- ggplot(data=data.frame(model_quantile,max_order),
                aes(x=max_order,y=model_quantile)) +
  geom_point(shape = 1, col = "#33666C") +
  geom_abline(intercept=0,slope=1,col="red") +
  theme_piss(16, 11) +
  labs(y = "Model quantile", x = "Empirical quantile") +
  ggtitle("QQ-plot")
# Same conclusion


gridExtra::grid.arrange(gg_qq,gg_pp, nrow = 1)



#############################
gev.diag(gev_tx) # confidence bands seem not too large in the return leel plot !


################## Return Level ######################
#####################################################

# 100-year return level ?
y_100 <- -(1 / log(1 - 1/100))
r_m100 <- gev_tx1$estimate[1] + (gev_tx1$estimate[2] / gev_tx1$estimate[3]) *
  (y_100^gev_tx1$estimate[3] - 1)
# see eq. in the text --> here we take the estimate
GEV.INV(1 - 1/100, gev_tx1$estimate[1], gev_tx1$estimate[2], gev_tx1$estimate[3])
# or directly by our Inverse function actually, at survival probability

extRemes::return.level(gev_tx2, return.period = c(2, 10, 100, 10000))
## Beware that these assume that data are stationary, while they are not (?) !!
# We remark clearly that these inferences on return levels do not take into account
# the climate warming, ie the fact that mean temperature increases slightly over time
# AND the fact that extremes are more frequent in a climate change [ref] (...) Indeed,
#one shall not expect a 100-year return level TX of 38deg reached only after 10k years....
rl_proflik <- extRemes::return.level(gev_tx2, return.period = c(2, 10, 100, 10000),
                                     method = "proflik", do.ci = F)
rl_proflik

extRemes::return.level(gev_tx2, return.period = 50,
                       do.ci = T, method = "proflik")
extRemes::return.level.fevd.mle(gev_tx2, return.period = 200,
                                do.ci = T, method = "proflik")

normal_rl_ci <- extRemes::return.level.fevd.mle(gev_tx2, return.period = c(2,5,10,100,1e3),
                                                method = "normal", do.ci = T)
normal_rl_ci

rl_proflik2 <- ci(gev_tx2,method="proflik",xrange=c(25, 38),verbose=T, return.period = 2)
rl_proflik2
rl_proflik5 <- ci(gev_tx2,method="proflik",xrange=c(30, 40),verbose=T, return.period = 5)
rl_proflik5
rl_proflik10 <- ci(gev_tx2,method="proflik",xrange=c(30, 40),verbose=T, return.period = 10)
rl_proflik10
rl_proflik100 <- ci(gev_tx2,method="proflik",xrange=c(35, 40),verbose = T)
rl_proflik100
rl_proflik1k <- ci(gev_tx2,method="proflik",xrange=c(35, 45),verbose=T, return.period = 1e3)
rl_proflik1k

# Estimates are the same as profiled likelihood...


## Standard errors of the estimates

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



#### Profile likelihood  ( For stationary model !!! )
par(mfrow = c(1,3))
gev.prof(gev_tx, 2, xlow = min(max_years$data)+2, xup = max(max_years$data)+5)
title(main= "2-year return level")
gev.prof(gev_tx, 10, xlow = min(max_years$data)+2, xup = max(max_years$data)+5)
title(main= "10-year return level")
gev.prof(gev_tx, 100, xlow = min(max_years$data)+2, xup = max(max_years$data)+5)
title(main= "100-year return level")
## Visualization is bad for high return levels. Let's rewrite the function...
# Here, we will use base plot because we don't want to loose too much time on this...

# This function allows to directly to compute the profile likelihood intervals plot for
# the return level, where y_1 and y_2 parameters allows to handle the y-coordinates
# in the plots.
"gev_prof_rl" <- function(mle, m, xlow, xup,
                          y_1 = 0, y_2 = 0, conf.ci = 0.95, n.int = 100) {
  #browser()
  z <<- mle  ;       p <<- 1/m
  v <- numeric(n.int)
  x <- seq(xlow, xup, length = n.int)
  sol <- c(z$mle[2], z$mle[3])

  for(i in 1:n.int) {
    xp <<- x[i]
    opt <- optim(sol, PissoortThesis::gev.proflik)
    sol <- opt$par ; v[i] <- - opt$value
  }
  plot(x,  v, type = "l", ylim = c(min(v) + y_1, max(v) + y_2),
       xlab = "Return Level", ylab ="Profile Log-likelihood")
  ma <-  - z$nllh
  #abline(h = ma, col = 4)
  abline(h = ma - 0.5 * qchisq(conf.ci, 1), col = 4)
}

## Make the multi-plot !
par(mfrow = c(1,3))
gev_prof_rl(gev_tx, 2, xlow = min(max_years$data)+2, xup = max(max_years$data)+5)
abline(v = rl_proflik2[1], col = 3)
abline(v = rl_proflik2[3], col = 3)
title(main= "2-year return level")
gev_prof_rl(gev_tx, 10, xlow = min(max_years$data)+2, xup = max(max_years$data)+5,
            y_1 = 150)
abline(v = rl_proflik10[1], col = 3)
abline(v = rl_proflik10[3], col = 3)
title(main= "10-year return level")
gev_prof_rl(gev_tx, 100, xlow = min(max_years$data)+2, xup = max(max_years$data)+5,
            y_1 = 450)
abline(v = rl_proflik100[1], col = 3)
abline(v = rl_proflik100[3], col = 3)
title(main= "100-year return level")

# Not accurate if too long-term prediction
gev.prof(gev_tx, 1000, xlow = min(max_years$data)+5, xup = max(max_years$data)+5)
gev.prof(gev_tx, 1e5, xlow = min(max_years$data)+5, xup = max(max_years$data)+10)




# Return Levels and empirical Quantiles

rl_df <- PissoortThesis::rl_piss(gev_tx$mle, gev_tx$cov, gev_tx$data)

gg_rl1 <- ggplot(rl_df) + coord_cartesian(xlim = c(0.1, 1000)) +
  scale_x_log10(breaks = c(0, 1,2, 5, 10,100, 1000), labels = c(0, 1, 2, 5, 10,100, 1000)) +
  labs(title = " Return Level plot", x = "Year (log scale)", y = "Quantile") +
  geom_point(aes(x = point, y = pdat), col = "#33666C", shape = 1) +
  geom_hline(yintercept = upp_endpoint, linetype = 2, col = "green3") +
  geom_hline(yintercept = max(max_years$data), linetype = 2, size = 0.25) +
  theme_piss()

intervals <- c( "Normal" = "red", "Prof.Likelihood" = "blue")
gg_rl <- gg_rl1 +
  ## Observations + Normal intervals
  geom_line(aes(x = y, y = q), col = "#33666C") +
  geom_line(aes( x = y, y = low, col = "Normal")) +
  geom_line(aes (x = y, y = upp, col = "Normal")) +
  # Profile likehood intervals (blue dots)
  geom_point(aes(x = 2, y = rl_proflik2[1], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 2, y = rl_proflik2[3], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 5, y = rl_proflik5[1], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 5, y = rl_proflik5[3], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 10, y = rl_proflik10[1], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 10, y = rl_proflik10[3], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 100, y = rl_proflik100[1], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 100, y = rl_proflik100[3], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 1e3, y = rl_proflik1k[1], col = "Prof.Likelihood"), size = 2) +
  geom_point(aes(x = 1e3, y = rl_proflik1k[3], col = "Prof.Likelihood"), size = 2) +
  scale_colour_manual(name = "Intervals", values = intervals,
                      guide = guide_legend(override.aes = list(
                          linetype = c("solid", "blank"),
                          shape = c(NA, 16),
                          size = c(1.4, 3)) ) ) +
  theme(legend.position = c(.89, .088)) +
  scale_y_continuous(breaks = c(30, 35, max(max_years$data), upp_endpoint, 40),
                     label = c("30", "35", as.character(round(max(max_years$data),1)),
                               as.character(round(upp_endpoint, 1) ), "40") )
  # geom_line(ymin = rl_proflik10[1], ymax = rl_proflik1000[1], xmin = 10, xmax = 100)
  # geom_area(aes(ymin = ))
gg_rl


## Density plots
x <- seq(min(max_years$data)-5, max(max_years$data)+5, length = length(max_years$data))
weib_fit_dens <- evd::dgev(x,loc = gev_tx$mle[1],
                           scale = gev_tx$mle[2], shape = gev_tx$mle[3])


density <- c( "empirical" = "black", "fitted" = "green3")
gg_ds <- ggplot(data.frame(x, weib_fit_dens, emp = max_years$data)) +
  stat_density(aes(x = emp, col = "empirical"),
               geom = "line", position = "identity") +
  #geom_density(aes (x = emp, col = "empirical")) +
  #geom_line(aes(x = x, y = emp, col = "empirical"))  +
  ggtitle("Empirical (black) vs fitted Weibull (green) density") +
  geom_line(aes(x = x, y = weib_fit_dens, col = "fitted"))  +
  coord_cartesian(xlim = c(25, 39)) +
  geom_vline(xintercept = min(max_years$data), linetype = 2) +
  geom_vline(xintercept = max(max_years$data), linetype = 2) +
  geom_vline(xintercept = upp_endpoint['loc'], linetype = 2, col = "green3") +
  theme_piss(17) + labs(x = "TX") +
  scale_colour_manual(name = "Density", values = density) +
  theme(legend.position = c(.915, .9))  +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(breaks = c( 25, min(max_years$data), 30, 35,
                                max(max_years$data), upp_endpoint),
                     label = c(25, as.character(round(min(max_years$data), 1)),
                               "30", "35",
                               as.character(round(max(max_years$data),1)),
                               as.character(round(upp_endpoint, 1) )) )
gg_ds
# Greeb is the fitted density while black is the empirical one.

gridExtra::grid.arrange(gg_rl, gg_ds, nrow = 1)

### TEST : Correct the visual 'problem' of the upper endpoints of kernel density
#===================================================================================
ggds_comp <- ggplot(max_years$df, aes (x = Max)) + geom_density() +
  geom_vline(xintercept = min(max_years$data), linetype = 2) +
  geom_vline(xintercept = max(max_years$data), linetype = 2) +
  coord_cartesian(xlim = c(25, 39))

gridExtra::grid.arrange(ggds_comp, gg_ds)


weib_fit_obs <- evd::rgev(116*800,loc = gev_tx$mle[1],
                          scale = gev_tx$mle[2], shape = gev_tx$mle[3])

dens <- density(max_years$df$Max, n = 2^7)
plot(dens)

ggplot() +  geom_line(aes(x = x, y = weib_fit_dens, col = "empirical"),
                      data = data.frame(x, weib_fit_dens))  +
  # stat_density(aes( col = "empirical"),
  #              geom = "line", position = "identity") +
  geom_density(data = max_years$df, aes (x = Max)) +
  ggtitle("Empirical (black) vs fitted Weibull (green) density") +
  # geom_density(aes(x = weib_fit_obs, col = "fitted"),
  #           data = data.frame(weib_fit_obs))  +
  coord_cartesian(xlim = c(25, 39)) + labs(x = "TX") +
  geom_vline(xintercept = min(max_years$data), linetype = 2) +
  geom_vline(xintercept = max(max_years$data), linetype = 2) +
  geom_vline(xintercept = upp_endpoint['loc'], linetype = 2, col = "green3") +
  theme_piss(17) +
  scale_colour_manual(name = "Density", values = density) +
  theme(legend.position = c(.915, .915))  +
  guides(colour = guide_legend(override.aes = list(size = 2)))


df <- data.frame(obs =  rep(max_years$data, 800), fit = weib_fit_obs )
df_melted <- reshape2::melt(df)
ggplot(df_melted, aes(x = value, fill = variable)) +
  geom_density(position = "stack", alpha = 0.6) +
  scale_y_continuous(name = "Density") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma")) +
  scale_fill_brewer(palette="Accent")
# ===============================================================================



##########################################################
###################    POT   #############################
##########################################################
library(POT)

max_all <- TXTN_closed$TX


### Choosing the threshold
##########################

## Mean Residual Life plot

# u <- seq(min(max_all),max(max_all),0.01)
# x <- c()
# for (i in 1:length(u)) {
#   threshold.excess <- max_all[max_all>u[i]]
#   x[i] <- mean(threshold.excess-u[i])
# }
# plot(x~u,type="l",main="MRL plot",ylab="mean excess", xlim = c(28,38), ylim = c(0,3))

## Several methods to do this plot, from packages or by hand ....
# mrlplot(max_all)
# mrl.plot(max_all)
mrl.plot(max_all, umin = 26, umax = max(max_all))  # Explore values
# Let's go for 31-32deg as threshold (?), as see ~some linearity in mrl plot

library(evmix)
mrlplot(max_all, tlim = c(20, max(max_all)))


## Selection based on stability of the estimates : TCplot
par(mfrow = (c(2,1)))
threshrange.plot(max_all, r = c(28,35) )
#tcplot(max_all, u.range = c(28, 35), nt = 50)

tcplot(max_all)# , tlim = c(25, max(max_all)-3),  )



### Do not consider the following pack as EVI<0 here. Just for testing
library(PissoortThesis)
MeanExcess(max_all, plot = T, k =1000)
genHill(max_all)
library(evir)
hill(max_all,start = 10,end = 10000)


## Dispersion Index Plot
# diplot(cbind(obs=TXTN_closed$t,time=TXTN_closed$Date,TXTN_closed$TX),
#        u.range = c(20,36) )
# Too slow... Remedy to this (were too much iterations in the loop)

data.di <- data.frame(time = seq(1:length(max_all)), obs = max_all)
events <- clust(data.di, u = 26, tim.cond = 20/365, clust.max = TRUE)
PissoortThesis::diplot_fast(events[,-3], u.range = c(25,35), nt = 1000 )

## L-moments plot
lmomplot(max_all, u.range = c(25, 35), identify = T)
# interpretation of this plot is really tedious.



#### Model Fit
###############

fit_pot <- fevd(max_all,threshold = 30, type="GP")   ## They use MLE by default
# pextRemes(fit_pot,c(25,30,35,36), lower.tail = FALSE)


gpd_mom <- fitgpd(max_all, 30, "moments") ;  gpd_mle <- fitgpd(max_all, 30, "mle")
gpd_pwmu <- fitgpd(max_all, 30, "pwmu") ;    gpd_pwmb <- fitgpd(max_all, 30, "pwmb")
gpd_pickands <- fitgpd(max_all, 30, "pickands")
gpd_med <- fitgpd(max_all, 30, "med", start = list(scale = 2, shape = 0.1))
gpd_mdpd <- fitgpd(max_all, 30, "mdpd") # mean power density divergence
gpd_mple <- fitgpd(max_all, 30, "mple")
gpd_ad2r <- fitgpd(max_all, 30, "mgf", stat = "AD2R") # maximum goodness-of-fit
#with modified Anderson Darling statistic
print(rbind(mom = gpd_mom$param, mle = gpd_mle$param,
            pwm.unbiased =  gpd_pwmu$param,
            pwm.biased = gpd_pwmb$param,
            pickands = gpd_pickands$param, median = gpd_med$param,
            mdpd =  gpd_mdpd$param, MpenalizedLE =  gpd_mple$param,
            ad2r = gpd_ad2r$param))

fit_pot1 <- gpd.fit(max_all, 30)
gpd.diag(fit_pot1)    #;   plot(fit_pot)
pextRemes(fit_pot,c(25,30,35,36), lower.tail = FALSE)


## (Fisher) Confidence intervals
# For the scale parameter
mle.ci <- gpd.fiscale(gpd_mle, conf = 0.9)
mom.ci <- gpd.fiscale(gpd_mom, conf = 0.9)
pwmu.ci <- gpd.fiscale(gpd_pwmu, conf = 0.9)
pwmb.ci <- gpd.fiscale(gpd_pwmb, conf = 0.9)
print(rbind(mle.ci,mom.ci,pwmu.ci,pwmb.ci))
# Shape parameter
shape_fi <- gpd.fishape(gpd_mle)
shape_pf <- gpd.pfshape(gpd_mle, range = c(-0.4, -0.1), conf = c(.95))
print(rbind(shape_fi,shape_pf))
# Profiled likehood interval is more narrow (...)


## Point Process
#################
gev_pp <- fevd(max_years$data, units = "deg C", threshold = 30,
               type = "PP", verbose = T)
plot(gev_pp)
## Beware : For non-stationary models, the density plot is not provided.
plot(gev_pp, "trace")
ci(gev_pp, type = "parameter")
threshrange.plot(max_all, r = c(27, 34), nint = 20, type = "PP")





# Not possible to do POT for whole days in year. max in winter will not be max etc ...#
#######################################################################################

# see distribution of TX wrt "arbitrary' threshold
above_25 <- TXTN_closed[TXTN_closed$TX>25,] #  ( change value of 25)

ggplot(data = as.data.frame(above_25), aes(x = as.factor(month), y = TX)) +
  geom_boxplot() + theme_piss()
# Only months from april to october have TX excessing 25degC


## Clusters ?
above_32 <- TXTN_closed[TXTN_closed$TX>32, ] #  ( change value of 32)
ggplot(above_32) + geom_point(aes(x = Date, y = TX))
# Indeed clusters of exceedance with POT as exceedances are not independent
# ---> stationary ?   One can also remark particular warm years, like 1976
# Declustering techniques ?



#################### Idead : Dividing analysis Winter-Summer, seasons, etc...  #######
######################################################################################


TXTN_closed_s <- TXTN_closed[TXTN_closed$month %in% c(4, 5, 6, 7, 8, 9), ]
TXTN_closed_w <- TXTN_closed[TXTN_closed$month %in% c(1, 2, 3, 10, 11, 12), ]

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


## Retrieve the pertaining datasets w and s, to enable compute TX/TN by "seasons" (w or s)
list_years_w <- split(TXTN_closed_w,TXTN_closed_w$year)
list_years_s <- split(TXTN_closed_s,TXTN_closed_s$year)



######    WINTER

max_w <- yearly.extrm(list_years_w)
# Replace the last value which is too low as TX for winter did not occur
max_w$df[116,]$Max <- mean(max_w$df$Max[110:115])

summary(lm_w <- lm(max_w$df$Max ~ max_w$df$Year))

ggplot(data = max_w$df, aes(x = Year,y = Max)) + geom_point() + theme_bw()
ggplot(data = max_w$df, aes(x = Max)) + geom_histogram() + theme_minimal()
ggplot(data = max_w$df, aes(x = Max)) + geom_density() + theme_bw()



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
## Winter

min_w <- yearly.extrm(list_years_w, Fun = min, tmp = 'TN')

summary(lm_w_min <- lm(min_w$df$Min ~ min_w$df$Year))

ggplot(data =min_w$df ,aes(x = Year, y = Min)) + geom_line() +
  geom_smooth(method = 'lm',formula = y~x) + theme_bw() +
  stat_smooth(method = "loess", se = F, col = 'red')

ggplot(data = min_w$df, aes(x = Year, y = Min)) + geom_point() + theme_bw()
ggplot(data = min_w$df, aes(x = Min)) + geom_histogram() + theme_minimal()
ggplot(data = min_w$df, aes(x = Min)) + geom_density() + theme_bw()

## Summer

min_s <- yearly.extrm(list_years_s, Fun = min, tmp = "TN")

summary(lm_s_min <- lm(min_s$df$Min ~ min_s$df$Year))

ggplot(data = min_s$df, aes(x = Year,y = Min)) + geom_point() + theme_bw()
ggplot(data = min_s$df, aes(x = Min)) + geom_histogram() + theme_minimal()
ggplot(data = min_s$df, aes(x = Min)) + geom_density() + theme_bw()



gmin.w <- ggplot(data =min_w$df ,aes(x = Year, y = Min)) +
  geom_line(size = 0.3) +
  geom_smooth(method = 'lm', formula = y~x, size = 0.5) +
  stat_smooth(method = "loess", se = F, col = 'red', size = 0.5) +
  ggtitle("TN for winter months [10-03]") +
  theme_piss(14, 12) +
  theme(legend.title = element_text(colour="#33666C",
                                    size=12, face="bold")) +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.y = element_text(size = 5))

gmin.s <- ggplot(data = min_s$df, aes(x = Year, y = Min)) +
  geom_line(size = 0.3) +
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
# Variance well more for minima in winter than in summer months.


### Or more convenient to work with NEGATED data for TN ?
neg_min_s <- min_s$df   ;   neg_min_s$negMin <- -neg_min_s$Min
neg_min_w <- min_w$df   ;   neg_min_w$negMin <- -neg_min_w$Min

gnegmin.w <- ggplot(data = neg_min_w ,aes(x = Year, y = negMin)) +
  geom_line() +
  geom_smooth(method = 'lm', formula = y~x) +
  theme_bw() +
  stat_smooth(method = "loess", se = F, col = 'red') +
  ggtitle("Min for winter months")
gnegmin.s <- ggplot(data = neg_min_s, aes(x = Year, y = negMin)) +
  geom_line() +
  geom_smooth(method = 'lm', formula = y~x) +
  theme_bw() +
  ggtitle("Min for summer months") +
  stat_smooth(method = "loess", se = F, col = 'red')

grid.arrange(gnegmin.w, gnegmin.s, nrow = 2)


#################
# As expected, the trend is (slightly) heavier for TX in warm months
#  and heavier for TN in cold month. Test if it is really significant
#############

summary(glm(TXTN_closed_ws$TX ~ max_s$df$Year * TXTN_closed_ws$period))
summary(lm(TXTN_closed_ws$TX~TXTN_closed_ws$period * TXTN_closed_ws$year))
summary(lm(TXTN_closed_ws$TX~TXTN_closed_ws$period))




##############################################################################
##################    POT taken on divided datasets   #######################
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

plot(above_thres_s[1:length(above_thres_s)-1],
     above_thres_s[2:length(above_thres_s)] )
# Very slight sign of dependence....


## Now look at the next code dealing with nonstationarity.  !!!!



#save.image("data1.Rdata") #  To load in the other scripts
