# load("/home/piss/Documents/Extreme/R resources/IRM/data1.Rdata")
# source('/home/piss/PissoortRepo/PissoortThesis/Scripts-R/Funs_introSplines.R')

load("/home/proto4426/Documents/master_Thesis/Extreme/R resources/IRM/data1.Rdata")
source('/home/proto4426/Documents/Thesis/PissoortThesis/Scripts-R/Funs_introSplines.R')
# Functions used that we will source. These are not our property so we will not include
# them in the package. These are mainly from Gavin Simpson (with some adjustments)


# Apply the created theme to all ggplots without having to specify it
theme_set(PissoortThesis::theme_piss())


library(tidyverse)
library(gridExtra)
library(PissoortThesis)
library(ggplot2)
library(mgcv)


## GAM(M) : Generalized Additive Models, Without any correlation structure from now.
gam1 <- gamm(Max ~ s(Year, k = 10), data = max_years$df)
# "this approach allows correlated errors to be dealt with via random effects or the
#corr structures available in nlme (using corr structures beyond the strictly additive case
#amounts to using a GEE approach to fitting).

summary(gam1$gam) # The smoother explains 20% of the variance of the data. 3 df

# look at autocorrelation in residuals:
acfpacf_plot(resid(gam1$lme, type = "normalized"), mfrow = c(1,2),
             main = "for the residuals of GAM1")
#  Hard to choose for AR or MA terms visually here...  should do further examinations (!)

gam2 <- gamm(Max ~ s(Year, k = 20), data = max_years$df,
             correlation = corARMA(form = ~ Year, p = 1, q = 0))
gam3 <- gamm(Max ~ s(Year, k = 20), data = max_years$df,
             correlation = corARMA(form = ~ Year, p = 0, q = 1))
gam4 <- gamm(Max ~ s(Year, k = 20), data = max_years$df,
             correlation = corARMA(form = ~ Year, p = 1, q = 1))
gam5 <- gamm(Max ~ s(Year, k = 20), data = max_years$df,
             correlation = corARMA(form = ~ Year, p = 2, q = 0))
gam6 <- gamm(Max ~ s(Year, k = 20), data = max_years$df,
             correlation = corARMA(form = ~ Year, p = 0, q = 2))

anovvv <- anova(gam1$lme, gam2$lme, gam3$lme, gam4$lme, gam5$lme, gam6$lme)
anovvv

library(stargazer)
latex <- cbind(df = anovvv$df, AIC = anovvv$AIC, BIC = anovvv$BIC)
rownames(latex) <- c("Indep", "AR(1)", "MA(1)", "ARMA(1,1)", "AR(2)", "MA(2)")
stargazer(latex)

anova(gam1$lme, gam3$lme)
anova(gam3$lme, gam4$lme)   # We would for for MA(1) : Simpler. Keep gam3 model


plot(gam3$gam, residuals = TRUE, pch = 19, cex = 0.75)
plot(gamm(Max ~ s(Year, k = 116), data = max_years$df,
     correlation = corARMA(form = ~ Year, p = 0, q = 1))$gam, residuals = TRUE, pch = 19, cex = 0.75)
summary(gam3$gam) # The smoother explains 20% of the variance of the data. 3 df

# Or with REML : can produce unbiased estimates of (co)variance parameters
#=> penalized spline model is expressed as a linear mixed model with wiggly bits of the spline treated as random effects,
gam3_REML <- gamm(Max ~ s(Year, k = 20), data = max_years$df,
             correlation = corARMA(form = ~ Year, p = 0, q = 1),
             method = "REML")
# Make the comparison !

# Diagnostic plots : (should redo it with ggplots)
with(max_years$df, tsDiagGamm(gam3, timevar = Year, observed = Max)) # See Fun.R



plot(Max ~ Year, data = max_years$df, type = "p", ylab = ylab)
pdat <- with(max_years$df, data.frame(Year = seq(min(Year),
                                                 max(Year),length = length(max_years$data))))
# Useful if lots of data and want to have a subset (pdat) of the data.
p1 <- predict(gam1$gam, newdata = pdat)
p2 <- predict(gam3$gam, newdata = pdat)
lines(p1 ~ Year, data = pdat, col = "red")
lines(p2 ~ Year, data = pdat, col = "blue")
legend("topleft",
       legend = c("Uncorrelated Errors","MA(1) Errors"),
       bty = "n", col = c("red","blue"), lty = 1)
# The MA(1) modelfor the residuals is questionnable... Model is not simpler.

## See if there there is a statistical difference in the different increasing patterns
# We estimate the first derivatives of the slope over a grid of points.
gam3.d <- Deriv(gam3$gam, n = length(max_years$data))  # See Fun.R
plot(gam3.d, sizer = TRUE, alpha = 0.01)# Shows first derivatives of the additive model with
#MA(1) errors, with 0.99 ci in grey area. Siggnificant increases are in blue.

# Surimpose it on the initial plot
plot(Max ~ Year, data = max_years$df, type = "p", ylab = 'Annual TX ' )
lines(p2 ~ Year, data = max_years$df)
CI <- confint(gam3.d, alpha = 0.01)
S <- signifD(p2, gam3.d$Year$deriv, CI$Year$upper, CI$Year$lower,
             eval = 0)
lines(S$incr ~ Year, data = max_years$df, lwd = 3, col = "blue")
lines(S$decr ~ Year, data = max_years$df, lwd = 3, col = "red")
# 2 periods of significant changes at 1% :


### 'improvements' of the method !

set.seed(10)
newd <- with(max_years$df, data.frame(Year = seq(min(Year), max(Year),
                                                 length.out = length(max_years$data))))
sims <- simulate(gam3, nsim = 10000, newdata = newd) # See Fun.R (??)

ci <- apply(sims, 1L, quantile, probs = c(0.025, 0.975))
newd <- transform(newd,
                  fitted = predict(gam1$gam, newdata = newd),
                  lower  = ci[1, ],
                  upper  = ci[2, ])

ggplot(max_years$df, aes(x = Year, y = Max)) + geom_point() +
  geom_ribbon(data = newd, aes(ymin = lower, ymax = upper, x = Year, y = fitted),
              alpha = 0.2, fill = "grey") +
  geom_line(data = newd, aes(y = fitted, x = Year))  + theme_piss()


# Get sense of the uncertainty in shapes of the simulatd trends : plots draws from posterior
set.seed(42)
S <- 50
sims2 <- setNames(data.frame(sims[, sample(1000, S)]), paste0("sim", seq_len(S)))
sims2 <- setNames(stack(sims2), c("TmaX", "Simulation"))
sims2 <- transform(sims2, Year = rep(newd$Year, S))
ggplot(sims2, aes(x = Year, y = TmaX, group = Simulation)) +
  geom_line(alpha = 0.3)
# illustrates the uncertainty in the estimates of the spline coefficients.
ggplot(sims2, aes(x = Year, y = TmaX, group = Simulation)) +
  geom_line(alpha = 0.5) + xlim(c(1975, 2016)) + ylim(c(30, 35))
# looks quite linear for this period !!


fd <- derivSimulCI(gam1, samples = 10000, n = 116) # See Fun.R

CI <- apply(fd[[1]]$simulations, 1, quantile, probs = c(0.025, 0.975))
sigD <- signifD(fd[["Year"]]$deriv, fd[["Year"]]$deriv, CI[2, ], CI[1, ],
                eval = 0)
newd <- transform(newd,
                  derivative = fd[["Year"]]$deriv[, 1], # computed first derivative
                  fdUpper = CI[2, ],                    # upper CI on first deriv
                  fdLower = CI[1, ],                    # lower CI on first deriv
                  increasing = sigD$incr,               # where is curve increasing?
                  decreasing = sigD$decr)               # ... or decreasing?

g_deriv_pointw <- ggplot(newd, aes(x = Year, y = derivative)) +
  geom_ribbon(aes(ymax = fdUpper, ymin = fdLower), alpha = 0.3, fill = "grey") +
  geom_line() + geom_hline(aes(yintercept = 0), col = "red") +
  geom_line(aes(y = increasing), size = 1.5) +
  geom_line(aes(y = decreasing), size = 1.5) +
  ylab(expression(hat(f) * "'" * (Year))) +
  labs(title = expression(paste(hat(f) * "'" * (Year), " of the fitted spline with .95 ", underline("pointwise"), " interval"))) +
  coord_cartesian(ylim = c(-0.09,0.22)) + theme_piss()
g_deriv_pointw


set.seed(123)
nsim <- 10000
pauseD <- derivSimulCI(gam1, samples = nsim,
                       newdata = data.frame(Year = seq(1975, 2016, by = 1)))

annualSlopes <- setNames(stack(setNames(data.frame(pauseD$Year$simulations),
                                     paste0("sim", seq_len(nsim)))),
                      c("Derivative", "Simulations"))
annualSlopes <- transform(annualSlopes, Year = rep(seq(1975, 2016, by = 1), each = nsim))

# Kernel density estimates of  first derivative of posterior sim from the fitted trend model for selected years
ggplot(annualSlopes, aes(x = Derivative, group = Year)) +
  geom_line(stat = "density", trim = TRUE) +
  geom_vline(xintercept = mean(annualSlopes$Derivative), col = "red") +
  theme_piss(theme = NULL, panel.background = element_rect(fill = NA),
             panel.grid.major = element_line(colour = "grey50"),
             panel.ontop = TRUE,
             axis.line = element_line(size = 3, colour = "grey80")) +
  facet_wrap(~ Year)

# smallest derivative for each year over all of the 10,000 posterior simulations
minD <- aggregate(Derivative ~ Year, data = annualSlopes, FUN = min)
ggplot(minD, aes(x = Year, y = Derivative)) +
  geom_point()
# Over 10k simulations, the isn't much which have derivatives <0...

library("viridis")
ggplot(annualSlopes, aes(x = Derivative, group = Year, colour = Year)) +
  geom_line(stat = "density", trim = TRUE) + scale_color_viridis(option = "magma") +
  theme(legend.position = "top", legend.key.width = unit(3, "cm")) + theme_piss()
# there’s little between-year shift in  slopes of the trends simulated from  posterior distribution of the model.
# similar to bayesian result of to Cahill et al. (2015)



# ==========================================================================
##  Identifying periods of sgnificant changes in the series with GAM

## fit an additive model with seasonal and trend smooth and an MA(1) process for the residuals;
#the code predicts from the model at 116 locations over the entire time series and generates
#a pointwise, approximate 95% confidence interval on the trend spline.
# Reuse above code (gam3, with MA(1) fitted on the errors)
n <- 116  # Try more !!!!!!!!! ( But it does not change anything for sure)
pdat <- with(max_years$df, data.frame(Year = seq(min(Year),max(Year),
                                                 length = n)))
p2 <- predict(gam3$gam, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])

df.res <- df.residual(gam3$gam)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))
head(pdat)  ; tail(pdat)

# method of finite differences to estimate derivative of a spline function in additive model


Term <- "Year"
gam3.d <- Deriv(gam3, n )  # see Fun.R  !!
# Why not 'tune', and inflate sample size to (perhaps) improve accuracy ??
gam3.dci <- confint.Deriv(gam3.d, term = Term) # OK if look only at the interval for
#single point in isolation but not if look at the entire spline => multiple comparisons issue
gam3.dsig <- signifD(pdat$p2, d = gam3.d[[Term]]$deriv,
                     gam3.dci[[Term]]$upper, gam3.dci[[Term]]$lower)
plot.Deriv(gam3.d)


ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(Temperature ~ (degree*C * ":" ~ centred))

plot(Max - mean(Max) ~ Year, data = max_years$df, type = "n",
     ylab = ylab, ylim = ylim)
points(Max - mean(Max) ~ Year, data = max_years$df,
        col = "lightgrey", pch = 16, cex = 0.7)
lines(p2 ~ Year, data = pdat)
lines(upper ~ Year, data = pdat, lty = "dashed")
lines(lower ~ Year, data = pdat, lty = "dashed")
lines(unlist(gam3.dsig$incr) ~ Year, data = pdat, col = "blue", lwd = 3)
lines(unlist(gam3.dsig$decr) ~ Year, data = pdat, col = "red", lwd = 3)
# >2deg increase in mean temperature since 1975, ...

df <- transform(max_years$df, TX_centered = Max - mean(Max), p2 = pdat$p2,
                upper = pdat$upper, lower = pdat$lower,
                increasing = unlist(gam3.dsig$incr), decreasing = unlist(gam3.dsig$decr))
ggplot(df, aes(x = Year)) +
  geom_point(aes(y = TX_centered), size = .9) +
  geom_line(aes(y = TX_centered), size = .55) +
  geom_line(aes(y = p2), size = 1, col = "green") +
  geom_line(aes(y = increasing), col = "darkgreen", size = 1.5) +
  geom_line(aes(y = decreasing), col = "darkgreen", size = 1.5) +
  geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.22, fill = "#33666C") +
  labs(title = "Centered Maxima with fitted GAM",
       y = expression(Centered ~ Max~(T~degree*C))) +
  theme_piss()



# =========================================================================


## Linear regression model... (Not very interesting actually : commented)


# ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
# m2 <- gamm(Temperature ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
#            data = cet, correlation = corARMA(form = ~ 1|Year, p = 2),
#            control = ctrl)
# summary(m2) ; m2$gam
# ## prediction data
# want <- seq(1, nrow(cet), length.out = 200)
# pdat <- with(cet,  data.frame(Time = Time[want], Date = Date[want],
#                               nMonth = nMonth[want]))
#
# set.seed(1)
# N <- 100
# dat <- data.frame(x = runif(N, min = 1, max = 20))
# dat <- transform(dat, y = 3 + (1.45 * x) + rnorm(N, mean = 2, sd = 3))
# ## sort dat on x to make things easier later
# dat <- dat[order(dat$x), ]
# summary(mod <- lm(y ~ x, data = dat) )
# (vc <- vcov(mod))
# coef(summary(mod))[, "Std. Error"] ;  sqrt(diag(vc))
#
# # Generate large number of values with mvnorm
# require("MASS")
# set.seed(10) ;   nsim <- 5000
# sim <- mvrnorm(nsim, mu = coef(mod), Sigma = vc) # Sigma here defines coVARIANCE matrix
# head(sim)
#
# kde <- kde2d(sim[,1], sim[,2], n = 75)
# plot(sim, pch = 19, cex = 0.5, col = "darkgrey")
# contour(kde$x, kde$y, kde$z, add = TRUE, col = "red",
#         lwd = 2, drawlabels = FALSE)
# # Large spread in the points illustrates the greater uncertainty in the
# # intercept term than in β^x. See also the negative correlation between these two.
# plot(y ~ x, data = dat)
# set.seed(42)
# take <- sample(nrow(sim), 50) ## take 50 simulations at random
# fits <- cbind(1, dat$x) %*% t(sim[take, ])
# matlines(dat$x, fits, col = "#A9A9A97D", lty = "solid", lwd = 2)
# abline(mod, col = "red", lwd = 1)
# matlines(dat$x, predict(mod, interval = "confidence")[,-1],
#          col = "red", lty = "dashed")
#


## Additive models

lp <- predict(gam3$gam, newdata = pdat, type = "lpmatrix")
coefs <- coef(gam3$gam) ;   vc <- vcov(gam3$gam)
# Generate few samples from the posterior
set.seed(35)   ;    n <- 25
sim <- mvrnorm(n, mu = coefs, Sigma = vc)
# Ignore terms relating to nMonth and intercept. Just want those pertaining to trend spline.
want <- grep("Year", colnames(lp))

fits <- lp[, want] %*% t(sim[, want])
dim(fits) ## n columns, 1 per simulation, 200 rows, 1 per evaln point.

# Draws
ylims <- range(fits)
plot(Max ~ Year, data = max_years$df, pch = 19, ylim = ylims, type = "n")
matlines(pdat$Year, fits, col = "black", lty = "solid")

# Posterior simulation for first derivatives of spline
fd <- derivSimulCI(gam3, samples = 10000, n = 116)

# take 250th and 9750th quantile for the ci if 10k resamples !
CI <- lapply(fd[1], function(x) apply(x$simulations, 1,
                                        quantile, probs = c(0.025, 0.975)))
plot.derivSimulCI(fd, sizer = TRUE)
# POINT-WISE confidence intervals !! see below, coverage is not 95%, so beware!

## Wrapping up

fit.fd <- fd[[1]]$deriv
set.seed(76)
take <- sample(nrow(fd[[1]]$simulations), 20) # use only 20 of the derivatives
plot(pdat$Year, fit.fd, type = "l", ylim = range(CI[[1]]), lwd = 2,
     main = "First derivative of trend spline from the series additive model.")
matlines(pdat$Year, t(CI[[1]]), lty = "dashed", col = "red")
matlines(pdat$Year, fd[[1]]$simulations[, take], lty = "solid",  col = "grey")
# => .95 pointwise c.i. with 20 of the samples surimposed.

## Beware that here, the c.i. were NOT simultaneous !! ==>




#  ===========================================================================
## SIMULTANEOUS intervals for smooths (revisited)


gam3.0 <- gam(Max ~ s(Year, k = 20), data = max_years$df, method = "REML")
# REMOVE correlation structure ! as it is not really significant ! Easier to handle.

summary(gam3.0) # ~3.5 (effective) df's. see ruppert(2003) p.81 for interpretation

plot(gam3.0, shade = TRUE, seWithMean = TRUE,   # see plot.gam()
     residuals = TRUE, pch = 16, cex = 0.8, rug = F,
     main = "fitted penalised spline with approximate 95% point-wise confidence interval")
# Shaded regions are confidence `bands` for smooths (not for parametric terms!).
# Notce that here, as we have same homegenous number of samples (1/year),
#confidence bands are likely to be the same everywhere (expcept at the edges).
# The confidence interval is a 95% Bayesian credible interval. For reasons
#This interval has a surprising frequentist interpretation as a 95%
#“across the function” interval; under repeated resampling from the population,
#95% of such confidence intervals will contain the true function.
# This is intuitive BUT does not reflect true uncertainty in the fitted function.
# Fewer than 95% of splines drawn from posterior of the fitted GAM would lie
#within the confidence interval shown in the plot above

# From Ruppert et al. (2003) `Semiparametric Regression´

# Now, rather than simulating from the posterior of the model as done before,
#this time we simulate from multivariate normal distribution with mean vector 0 and
#covariance matrix Vb, the Bayesian covariance matrix of the fitted model

# to generate random values from a multivariate normal
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig) ;   m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

Vb <- vcov(gam3.0) # Bayesian covariance matrix of the model coefficients.
#it's conditional upon the smoothing parameter(s). add unconditional=T
#to adjust for the smoothing parameters being estimated rather than known values,

grid <- 116 # Defines grid of values for which we will generate predictions.
newd <- with(max_years$df,
             data.frame(Year = seq(min(Year), max(Year), length = grid)))
pred <- predict(gam3.0, newd, se.fit = TRUE)
se.fit <- pred$se.fit


set.seed(42)
N <- 10000
# We want N draws from [\hat{beta}-beta, \hat{u}-u] which is ~multi N(0,Vb)
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

# Then compute \hat{f}(x)-f(x) which is  C_g%*%[\hat{beta}-beta, \hat{u}-u]
Cg <- predict(gam3.0, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

# Find absolute values of the standardized deviations from the true model
absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))

# Max of the absolute standardized dev. at the grid of x values for each simul
masd <- apply(absDev, 2L, max)

# Find the crit value used to scale standard errors to yield the simultaneous interval
crit <- quantile(masd, prob = 0.95, type = 8) # = 2.9
# compared with pointwise, it is now 2.9/1.96 ~ 1.5 times wider !


# Now, compute and show the simultaneous confidence interval !
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (qnorm(.975) * se.fit),
                  lwrP = fit - (qnorm(.975) * se.fit),
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))
ggplot(pred, aes(x = Year)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "red") +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2, fill = "red") +
  ggtitle("Comparison of point-wise and simultaneous 95% conf intervals for fitted GAM")+
  theme_piss()


# Look at the coverage properties of this interval ! If we do that for the previous
# (pointwise) intervals computed, we will see an incorrect coverage !

sims <- rmvn(N, mu = coef(gam3.0), sig = Vb)
fits <- Cg %*% t(sims) # contains N draws from the posterior

# First, we choose nrnd samples to represent it
nrnd <- 50   ;    rnd <- sample(N, nrnd)
stackFits <- stack(as.data.frame(fits[, rnd]))
stackFits <- transform(stackFits, Year = rep(newd$Year, length(rnd)))

interval <- c("pointwise" =  "yellow", "simultaneous" = "darkred")


ggplot(pred, aes(x = Year, y = fit)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = "simultaneous"), alpha = 0.4) +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP, fill = "pointwise"), alpha = 0.4) +
  geom_path(lwd = 2) +
  geom_path(data = stackFits, mapping = aes(y = values, x = Year, group = ind),
            alpha = 0.5, colour = "grey20") +
  labs(y = expression( Max~(T~degree*C)), x = "Year",
       title = "Point-wise & Simultaneous 95% conf. intervals for fitted GAM",
       subtitle = sprintf("Each line is one of %i draws from the posterior distribution of the model", nrnd)) +
  annotate(geom = "text", label = paste("coverages", " are : \n",
                                        round(pointw_cov, 5),
                                        " for pointwise \n", "   ",
                                        round(simult_cov, 5),
                                        " for simultaneous"),
           x = 1915,
           y = 33.4, col = "#33666C" , size = 4) +
  scale_fill_manual(name = "Interval", values = interval) +
  theme_piss(size_p = 15, plot.subtitle = text(12, hjust = 0.5, colour = "#33666C"), legend.position = c(.888, .152)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
# Many lines lie outside the pointwise interval !


## Compute the coverages

'inCI' <- function(x, upr, lwr) {
  # =T if all evaluation points g lie within the interval and =F otherwise.
  all(x >= lwr & x <= upr)
}

fitsInPCI <- apply(fits, 2L, inCI, upr = pred$uprP, lwr = pred$lwrP)
fitsInSCI <- apply(fits, 2L, inCI, upr = pred$uprS, lwr = pred$lwrS)

pointw_cov <- sum(fitsInPCI) / length(fitsInPCI)  # Point-wise
simult_cov <- sum(fitsInSCI) / length(fitsInSCI)  # Simultaneous
# As expected, poitwise is .62<0.95 and simultaneous is 0.955 !

ddddff <- cbind.data.frame(N = 10000, "Pointwise" = pointw_cov, "Simultaneous" = simult_cov)
rownames(ddddff) <- ("Coverage at 95%")
stargazer(ddddff)

oldCI <- apply(fits, 1L, quantile, probs = c(0.025, 0.975))
pred <- transform(pred, lwrOld = oldCI[1, ], uprOld = oldCI[2, ])
fitsInOldCI <- apply(fits, 2L, inCI, upr = pred$uprOld, lwr = pred$lwrOld)

sum(fitsInOldCI) / length(fitsInOldCI)
# With the CI computed above, it is lessthan 60% coverage !!

# We have computed actual simultaneous interval but only for the fitted smoother


ggplot(pred, aes(x = Year, y = fit)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = "simultaneous"), alpha = 0.4) +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP, fill = "pointwise"), alpha = 0.4) +
  geom_path(lwd = 2) +
  geom_path(data = stackFits, mapping = aes(y = values, x = Year, group = ind),
            alpha = 0.5, colour = "grey20") +
  labs(y = expression( Max~(T~degree*C)), x = "Year",
       title = "Point-wise & Simultaneous 95% conf. intervals for fitted GAM",
       subtitle = sprintf("Each line is one of %i draws from the posterior distribution of the model", nrnd)) +
  scale_fill_manual(name = "Interval", values = interval) +
  theme_piss(size_p = 15, plot.subtitle = text(12, hjust = 0.5, colour = "#33666C"), legend.position = c(.888, .152)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))


#  ===========================================================================

# Now,we will do similar thing but applied to the derivatives of a spline.
# We keep the gam3.0 model defined above

#devtools::install_github("gavinsimpson/tsgam")
library(tsgam) # see it


## parameters for testing
UNCONDITIONAL <- FALSE # conditional on estimating smooth params
#correction uncertainty due to smoothing parameters being estimated (not done here)
N <- 10000             # number of posterior draws
n <- 500               # number of newdata values. This number is not very imoportant here
EPS <- 1e-07           # finite difference

## where are we going to predict at?
newd <- with(max_years$df,
             data.frame(Year = seq(min(Year), max(Year), length = n)))
# See tsgam package
fd <- fderiv(gam3.0, newdata = newd, eps = EPS, unconditional = UNCONDITIONAL)
# computes the first derivative of any splines in the supplied gam(m)
str(fd, max = 1)

set.seed(42)
head( sint <- confint(fd, type = "simultaneous", nsim = N) )

g_deriv_simult <-
  ggplot(cbind(sint, Year = newd$Year), aes(x = Year, y = est)) +
  geom_hline(aes(yintercept = 0), col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() + coord_cartesian(ylim = c(-0.09,0.22)) +
  ylab(expression(hat(f) * "'" * (Year))) +
  labs(title = expression(paste(hat(f) * "'" * (Year), " of the fitted spline with .95 ", underline("simultaneous"), " interval"))) +
  theme_piss()
g_deriv_simult
# Much wider than those obtained above !! See it
gridExtra::grid.arrange(g_deriv_pointw, g_deriv_simult, nrow = 1)


# all we need to do to apply it to derivatives is to make the assumption that
#estimate of the first derivative is unbiased and hence we can proceed as we did previously
tsgam:::simultaneous

Vb <- vcov(gam3.0, unconditional = UNCONDITIONAL)
set.seed(24)
sims <- MASS::mvrnorm(N, mu = coef(gam3.0), Sigma = Vb)
X0 <- predict(gam3.0, newd, type = "lpmatrix")
newd <- newd + EPS
X1 <- predict(gam3.0, newd, type = "lpmatrix")
Xp <- (X1 - X0) / EPS
derivs <- Xp %*% t(sims)

set.seed(2)
matplot(derivs[, sample(N, 50)], type = "l", lty = "solid")
# Notice also the large uncertainty...


## Coverage properties : For the derivative!!  (Reuse inCI function defined above)

fitsInCI <- with(sint,
                 apply(derivs, 2L, inCI, upr = upper, lwr = lower))
sum(fitsInCI) / length(fitsInCI)  # Yeaahhhh  0.95 !!



