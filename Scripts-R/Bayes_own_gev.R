setwd('/home/piss/Documents/Extreme/R resources/IRM')
library(evd)
library(mvtnorm)
library(KernSmooth)
library(coda)
library(pander)
load("data1.Rdata")

library(PissoortThesis)

############   MEtropolis-Hastlings    ###############
#####################################################


# Optimize Posterior Density Function to find starting values
fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
param <- c(mean(max_years$df$Max),log(sd(max_years$df$Max)), 0.1 )
# opt <- optim(param, fn, data = max_years$data,
#              method="BFGS", hessian = TRUE)
opt <- nlm(fn, param, data = max_years$data,
           hessian=T, iterlim = 1e5)
start <- opt$estimate
Sig <- solve(opt$hessian)
ev <- eigen( (2.4/sqrt(2))^2 * Sig)
varmat <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)

set.seed(100)
mh.mcmc1 <- MH_mcmc.own(start, varmat %*% c(.1,.3,.4))
mh.mcmc1$mean.acc_rates

chains.plotOwn(mh.mcmc1$out.chain)





##########  GIBBS sampler  #####################
###############################################

#prop_var <- sqrt( (2.4/sqrt(1))^2 * solve(opt$hessian) )

set.seed(100)
iter <- 2000
gibb1 <- gibbs_mcmc.own(start, iter = iter) # Same starting point as MH
gibb1$mean.acc_rates

# Do not forget Burn in period (We will make it inside function in  following)
burn <- iter/4  # Tune value

gibb1$out.chain <- gibb1$out.chain[-(1:burn),]


chains.plotOwn(gibb1$out.chain)


param_gibbs <- apply(gibb1$out.chain[,1:3], 2, mean)
param_gibbs["logsig"] <- exp(param_gibbs["logsig"] )

frame <- data.frame(Bayesian = param_gibbs, 'Frequentist(mle)' = gev_fit$mle)
row.names(frame) = c("$\\mu \\ $", "$\\sigma \\quad$", "$\\xi \\quad$")
knitr::kable(frame, align = 'l')



#####################################################################
################### Gibbs Sampler with Nonstationarity ##############

data <- max_years$data


fn <- function(par, data) -log_post1(par[1], par[2], par[3],
                                        par[4], data)
param <- c(mean(max_years$df$Max), 0, log(sd(max_years$df$Max)), -0.1 )
opt <- optim(param, fn, data = max_years$data,
             method = "BFGS", hessian = T)
# opt <- nlm(fn, param, data = max_years$data,
#     hessian=T, iterlim = 1e5)
opt


# Starting Values
set.seed(100)
start <- list() ; k <- 1
while(k < 5) { # starting value is randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt$par, 50 * solve(opt$hessian)))
  svlp <- log_post1(sv[1], sv[2], sv[3], sv[4], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start[[k]] <- sv
    k <- k + 1
  }
}
knitr::kable(matrix(unlist(start), ncol = 4, byrow = T,
                    dimnames = list(c("start[[1]]", "start[[2]]",
                                      "start[[3]]", "start[[4]]"),
                                    c("$\\mu \\ $","$\\mu_{trend}$",
                                      "$\\logsig$", "$\\xi$"))), align = "c")


# k chains with k different starting values
set.seed(100)
gibbs.trend <- gibbs.trend.own(start, propsd = c(.5, 1.9, .15, .12),
                               iter = 1000)
colMeans(do.call(rbind, gibbs.trend$mean_acc.rates))


param.chain <- gibbs.trend$out.chain[ ,1:4]



### Plot of the chains
chains.plotOwn(gibbs.trend$out.chain )
ggplot(gibbs.trend$out.chain) + geom_line(aes(x = iter, y = mu1)) +
  theme_piss(16,14) + labs(ylab = "mu1", xlab = "iter" )


## TracePlots
chain.mix <- cbind.data.frame(gibbs.trend$out.chain,
                              iter.chain = rep(1:500, 4))
mixchains.Own(chain.mix)
ggplot(chain.mix, aes(x = iter.chain, y = mu1, col = as.factor(chain.nbr))) +
  geom_line() + theme_piss(18,16, theme_classic()) +
  scale_colour_brewer(name = "chain nr", palette = "Set1") +
  guides(colour = guide_legend(override.aes = list(size= 1.2)))

# Bayes Plots
array.post <- array(unlist(gibbs.trend$out.ind), dim = c(4000/4+1, 4, 4 ),
                    dimnames = list(iterations = NULL,
                                    parameters = c("mu", "mu1", "logsig", "xi"),
                                    chains = c("chain:1", "chain:2", "chain:3",
                                               "chain:4") ))

color_scheme_set("mix-blue-red")
mcmc_trace(array.post, pars = c("mu", "mu1"),
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_trace_highlight(array.post, pars = "mu", highlight = 3)


# Function to create mcmc.lists, useful for diagnostics on chains.
mc.listDiag <- function (list, subset = c("mu", "mu1", "logsig", "xi")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset]),
            mcmc(list[[4]][, subset])
            )
}
## Gelman Coda Diagnostics
gelman.diag(mc.listDiag(gibbs.trend$out.ind), autoburnin=F)
gelman.plot(mc.listDiag(gibbs.trend$out.ind), autoburnin=F)



# Markov Chain Correlations
autocorr(mcmc(param.chain ))
autocorr.diag(mcmc(param.chain ))
autocorr.plot(mcmc(param.chain ))

crosscorr(mcmc(param.chain))
crosscorr.plot(mcmc(param.chain ))


## Geweke
geweke <- geweke.diag(mcmc(param.chain))
2*dnorm(geweke$z)
geweke.plot(mcmc(param.chain), nbins = 20)



# Raftery Coda Diagnostics
raftery.diag(mc.listDiag(gibbs.trend$out.ind)[[1]], q=0.05, r=0.02, s=0.95)
raftery.diag(mc.listDiag(gibbs.trend$out.ind)[[2]], q=0.05, r=0.02, s=0.95)
raftery.diag(mc.listDiag(gibbs.trend$out.ind)[[3]], q=0.05, r=0.02, s=0.95)
raftery.diag(mc.listDiag(gibbs.trend$out.ind)[[4]], q=0.05, r=0.02, s=0.95)


param.chain[, "sigma"] <- exp(param.chain[,"logsig"])
# Summary And Parameter Table
tab_quantiles <- as.data.frame(summary(mcmc(param.chain))$quantiles)
# colnames(tab_quantiles) <- c("$\\boldsymbol{q_{0.025}}$","$\\boldsymbol{q_{0.25}}$",
#                     "Median","$\\boldsymbol{q_{0.75}}$",
#                     "$\\boldsymbol{q_{0.975}}$")
# rownames(tab_quantiles) <- c("$\\mu \\ $","$\\mu_1 \\quad$", "$\\log\\sigma \\quad$",
#                     "$\\xi \\quad$", "$\\sigma$")

pander(tab_quantiles,split.table = Inf)


mean.post <- apply(param.chain, 2, mean)



## Comparisons with Frequentist's results (GEV)
par_gibbs_trend <- apply(gibbs.trend$out.chain[,c("mu", "mu1", "logsig", "xi")],
                         2, mean) # Average over the (3) generated chains
par_gibbs_trend["sigma"] <- exp(par_gibbs_trend["logsig"] )
names(par_gibbs_trend["logsig"]) <- "sigma"
frame <- data.frame(Bayesian = par_gibbs_trend,
                    'Frequentist(mle)' = gev_nonstatio$mle)
row.names(frame) = c("$\\mu \\ $","$\\boxed{\\mu_1} \\quad$",
                     "$\\log\\sigma \\quad$", "$\\xi \\quad$", "$\\sigma$")
knitr::kable(frame, align = 'l')

## Credible intervals (quantile-based) see above quantiles ! Densities are quite
# symetric, so it is not a bad idea to use this method.

library("bayesplot")
color_scheme_set("blue")
mcmc_areas(
  posterior,
  pars = c("cyl", "drat", "am", "sigma"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)
mcmc_intervals(param.chain, pars = c("logsig", "xi"))


## Densities of the parameters with their quantile-based 0.95 intervals
color_scheme_set("brightblue")
g1 <- mcmc_dens(param.chain, pars = c("mu")) +
  geom_vline(xintercept = tab_quantiles['mu', "2.5%"], col = "red") +
  geom_vline(xintercept = tab_quantiles['mu', "97.5%"], col = "red")
g2 <- mcmc_dens(param.chain, pars = c("sigma")) +
  geom_vline(xintercept = tab_quantiles['sigma', "2.5%"], col = "red") +
  geom_vline(xintercept = tab_quantiles['sigma', "97.5%"], col = "red")
g3 <- mcmc_dens(param.chain, pars = c("xi"))+
  geom_vline(xintercept = tab_quantiles['xi', '2.5%'], col = "red") +
  geom_vline(xintercept = tab_quantiles['xi', "97.5%"], col = "red")
g4 <- mcmc_dens(param.chain, pars = c("mu1")) +
  geom_vline(xintercept = tab_quantiles['mu1', '2.5%'], col = "red") +
  geom_vline(xintercept = tab_quantiles['mu1', "97.5%"], col = "red")
grid.arrange( g1, g2, g3, g4, nrow = 2 )


## HPD intervals
library(HDInterval)
hpd_mu <- hdi(param.chain$mu)
hpd_sigma <- hdi(param.chain$sigma)
hpd_xi <- hdi(param.chain$xi)
hpd_mu1 <- hdi(param.chain$mu1)

### Comparisons of all the ci :







## Problem handling of mu_trend
t_bayes <- round(( min(max_years$df$Year):max(max_years$df$Year) -
                     mean(max_years$df$Year) ) / length(max_years$data),4)
t_freq <- seq(1, length(max_years$data),1)

mut <- par_gibbs_trend["mu"] + par_gibbs_trend["mu1"] * t_bayes
mu_freq <- gev_nonstatio$mle[1] + gev_nonstatio$mle[2] * t_freq

pander(t(matrix(c(mut[1:4], mut[(length(mut)-3):length(mut)],
                  mu_freq[1:4], mu_freq[(length(mu_freq)-3):length(mu_freq)]),
                ncol = 2, dimnames = list(c(rep("First values", 4),
                                            rep("Last values", 4)),
                                          c("Frequentist",  "Bayesian")))))


# http://www.biostat.umn.edu/~ph7440/


## Posterior Predictive Distribution
tt <- ( min(max_years$df$Year):(max(max_years$df$Year) ) -
          mean(max_years$df$Year) ) / length(max_years$data)

repl <- matrix(NA, nrow(gibbs.trend$out.chain), length(tt))
for(t in 1:nrow(repl)) {
  mu <- gibbs.trend$out.chain[t,1] + gibbs.trend$out.chain[t,2] * tt
  repl[t,] <- evd::rgev(length(tt), loc = mu,
                        scale = gibbs.trend$out.chain[t,3],
                        shape = gibbs.trend$out.chain[t,4])
}
post.pred <- apply(repl, 2, function(x) quantile(x, probs = c(0.05,0.5,0.95)))


df.postpred <- data.frame(data = max_years$data, q05 = post.pred["5%",],
                          q50 = post.pred["50%",], q95 = post.pred["95%",],
                          year = seq(1901:2016))
ggplot(df.postpred) + geom_point(aes(x = year, y = data)) +
  geom_line(aes(x = year, y = q05)) + geom_line(aes(x = year, y = q50)) +
  geom_line(aes(x = year, y = q95)) +
  ggtitle("Original data with PPD quantiles 5, 50 and 95%") + theme_piss()


## And for predictions ? 10 years here

tt2 <- ( min(max_years$df$Year):(max(max_years$df$Year) +10 ) -
           mean(max_years$df$Year) ) / length(max_years$data)

repl2 <- matrix(NA, nrow(gibbs.trend$out.chain), length(tt2))
for(t in 1:nrow(repl2)) {
  mu <- gibbs.trend$out.chain[t,1] + gibbs.trend$out.chain[t,2] * tt2
  repl2[t,] <- evd::rgev(length(tt2), loc = mu, scale = gibbs.trend$out.chain[t,3],
                         #shape = unname(par_gibbs_trend["xi"]))
                         shape = gibbs.trend$out.chain[t,4])
}
post.pred2 <- apply(repl2, 2, function(x) quantile(x, probs = c(0.05,0.5,0.95)))

df.postpred2 <- data.frame(org.data = c(max_years$data,
                                        repl2[sample(10, 1:nrow(repl2)),117:126] ),
                           q05 = post.pred2["5%",], q50 = post.pred2["50%",],
                           q95 = post.pred2["95%",], year = 1901:2026,
                           'data' = c(rep('original', 116), rep('new', 10)))

ggplot(df.postpred2, aes(col = data)) + geom_point(aes(x = year, y = org.data)) +
  geom_line(aes(x = year, y = q05)) + geom_line(aes(x = year, y = q50)) +
  geom_line(aes(x = year, y = q95)) +
  ggtitle("PPD quantiles 5, 50 and 95% and 10y predictions") + theme_piss()


hpd_pred <- hdi(repl2)

# Densities associated with the PPD, with mean(red) and
g1 <- ggplot(data.frame(repl2)) + stat_density(aes(x = X1), geom = "line") +
  labs(x = " TX for year 1901") + theme_piss() +
  geom_vline(xintercept = mean(repl2[,1]), col = "blue") +
  geom_vline(xintercept = post.pred2["5%", 1], col = "blue") +
  geom_vline(xintercept = post.pred2["95%", 1], col = "blue") +
  coord_cartesian(ylim = c(0.01,.55), xlim = c(27.5,36)) +
  geom_vline(xintercept = hpd_pred['lower', 1], col = "green") +
  geom_vline(xintercept = hpd_pred['upper', 1], col = "green")
g2 <- ggplot(data.frame(repl2)) + stat_density(aes(x = X50), geom = "line") +
  labs(x = " TX for year 1950") + theme_piss() +
  geom_vline(xintercept = mean(repl2[,50]), col = "blue") +
  geom_vline(xintercept = post.pred2["5%", 50], col = "blue") +
  geom_vline(xintercept = post.pred2["95%", 50], col = "blue") +
  coord_cartesian(ylim = c(0.01,.55), xlim = c(27.5,36)) +
  geom_vline(xintercept = hpd_pred['lower', 50], col = "green") +
  geom_vline(xintercept = hpd_pred['upper', 50], col = "green")
g3 <- ggplot(data.frame(repl2)) + stat_density(aes(x = X116), geom = "line") +
  labs(x = " TX for year 2016") + theme_piss() +
  geom_vline(xintercept = mean(repl2[,116]), col = "blue") +
  geom_vline(xintercept = post.pred2["5%", 116], col = "blue") +
  geom_vline(xintercept = post.pred2["95%", 116], col = "blue") +
  coord_cartesian(ylim = c(0.01,.55), xlim = c(27.5,36)) +
  geom_vline(xintercept = hpd_pred['lower', 116], col = "green") +
  geom_vline(xintercept = hpd_pred['upper', 116], col = "green")
g4 <- ggplot(data.frame(repl2)) + stat_density(aes(x = X126), geom = "line") +
  labs(x = " TX for year 2026 (?)") + theme_piss() +
  geom_vline(xintercept = mean(repl2[,126]), col = "blue") +
  geom_vline(xintercept = post.pred2["5%", 126], col = "blue") +
  geom_vline(xintercept = post.pred2["95%", 126], col = "blue") +
  coord_cartesian(ylim = c(0.01,.55), xlim = c(27.5,35)) +
  geom_vline(xintercept = hpd_pred['lower', 126], col = "green") +
  geom_vline(xintercept = hpd_pred['upper', 126], col = "green")
grid.arrange(g1, g2, g3, g4, nrow = 1)



## Predictive accuracy criterion

ic_vals <- gibbs.trend$dic.vals
# DIC Values
dic( mc.listDiag(gibbs.trend$out.ind)[[1]], ic_vals[[1]] )
dic( mc.listDiag(gibbs.trend$out.ind)[[2]], ic_vals[[2]] )
dic( mc.listDiag(gibbs.trend$out.ind)[[3]], ic_vals[[3]] )
dic( mc.listDiag(gibbs.trend$out.ind)[[4]], ic_vals[[4]] )
#handles uncertainty in inferences within each model, and it does
#not depend on aspects of the models that don’t affect inferences
#within each model. take a lot of simulations to calculate it precisely.
# despite the lack of a clear theoretical foundation.
# DIC is shown to be an approximation to a penalized loss function based
# on the deviance, with a penalty derived from a cross-validation argument.
# This approximation is valid only when the effective number of parameters
# in the model is much smaller than the number of independent observations

# WAIC Values
waic( ic_vals[[1]] )
waic( ic_vals[[2]])
waic( ic_vals[[3]])
waic( ic_vals[[4]] )


## Cross-validation



## Return Levels

# Stationary
gibb1$out.chain[, "sigma"] <- exp(gibb1$out.chain[, "logsig"])
rl.post_gg(gibb1$out.chain[, c("mu", "sigma", "xi")], method = "gev")
rl.pred_gg(gibb1$out.chain[, c("mu", "sigma", "xi")], method = "gev",
           qlim = c(30, 40), period = c(1, 5, 15))



# Nonstationary (linear trend)
library(evdbayes)
gibbs.trend$out.chain[,"sigma"] <- exp(gibbs.trend$out.chain[,"logsig"])
rl.data <-gibbs.trend$out.chain[, c("mu", "mu1", "sigma", "xi")]
rl.pred(rl.data, qlim = c(30, 40))

rl_bayes_trend <- function(data_year, params, t = 10, m = 10 ){
    y_m <- -(1 / log(1 - 1/m))
    t <- seq(max(data_year), max(data_year) + t, 1)
    rl_m <- (params[1] + params[2] *
               (t-max(max_years$df$Year))) +
      (params[3] / params[4]) *
      (y_m^params[4] - 1)
    g <- ggplot(data.frame(r.lvels = rl_m, years = t)) +
      geom_point(aes(x = years, y = r.lvels))
     # print(g)
    return(rl_m)
}
par_gibbs_trend[3] <- exp(par_gibbs_trend["logsig"])
rl_bayes_trend(max_years$data,
               unname(par_gibbs_trend[c("mu", "mu1", "sigma", "xi")]))





#### See end of bayesian01 for predictive dist, quantiles, and other...
