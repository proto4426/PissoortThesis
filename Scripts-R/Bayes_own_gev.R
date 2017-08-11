#setwd('/home/proto4426/Documents/Extreme/R resources/IRM')
#load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1.Rdata")
load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1_bayes.Rdata")


# Bayesian Analysis constrcuted with functions in package PissoortThesis
#######################################################################
library(evd)
library(mvtnorm)
library(KernSmooth)
library(coda)
library(pander)
library(gridExtra)
library(tidyverse)
library(ggjoy)
library(viridis)
library(grid)

library(PissoortThesis)



## As we aim at building models sequentially,  we will begin by the Gumbel model
## It is models the lowest number of degrees of freedom.



## MEtropolis-Hastlings for the stionary Gumbel model ###
#########################################################


# Optimize Posterior Density Function to find starting values
fn <- function(par, data) -log_post_gumb(par[1], par[2], data)
param <- c(mean(max_years$df$Max), log(sd(max_years$df$Max)))
opt <- nlm(fn, param, data = max_years$data,
           hessian=T, iterlim = 1e5)
start <- opt$estimate  ;     names(start) <- c("mu", "logsig")
Sig <- solve(opt$hessian)
ev <- eigen( (2.4/sqrt(2))^2 * Sig)
varmat <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)

set.seed(101)
iter <- 2e3
mh.mcmc_gumb <- PissoortThesis::MH_mcmc.own(start, varmat %*% c(.9, 1.1),
                                           iter = iter, burnin = iter/4)
mh.mcmc_gumb$mean.acc_rates



#########   MEtropolis-Hastlings for the stationary GEV model #######
#####################################################################


# Optimize Posterior Density Function to find starting values
fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
param <- c(mean(max_years$df$Max),log(sd(max_years$df$Max)), 0.1 )
# opt <- optim(param, fn, data = max_years$data,
#              method="BFGS", hessian = TRUE)
opt <- nlm(fn, param, data = max_years$data,
           hessian=T, iterlim = 1e5)
start <- opt$estimate  ;     names(start) <- c("mu", "logsig", "xi")
Sig <- solve(opt$hessian)
ev <- eigen( (2.4/sqrt(2))^2 * Sig)
varmat <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)

set.seed(101)
iter <- 2e3
mh.mcmc1 <- PissoortThesis::MH_mcmc.own(start, varmat %*% c(.5,.6,.85),
                                        iter = iter, burnin = iter/4)
mh.mcmc1$mean.acc_rates

param_mean_mh <- apply(mh.mcmc1$out.chain[,1:3], 2, mean)
chainsPlot_mh <- chains.plotOwn(mh.mcmc1$out.chain, post.mean.green = param_mean_mh,
                                title = "Using Metropolis-Hastings Algorithm")

## Burn in 1/4 of values  :
mh.mcmc.out <- mh.mcmc1$out.chain[-(1:iter/4), ]

# Effective sample sizes :
effectiveSize(mcmc(mh.mcmc.out[,1:3]))




##########  GIBBS sampler for the stationary GEV model  #######
##############################################################

#prop_var <- sqrt( (2.4/sqrt(1))^2 * solve(opt$hessian) )

set.seed(100)
iter <- 2e3
# Same starting point as for MH. We also keep only one chain
# i.e. no different starting values so far.
gibbs_statio <- PissoortThesis::gibbs_mcmc.own(start, iter = iter, nbr.chain = 1,
                                               propsd = c(.42, .12, .12),
                                               burnin = iter/4)
gibbs_statio$mean_acc.rates


chainsPlot_gibb <- chains.plotOwn(gibbs_statio$out.chain,
                                  title = "Using Gibbs Sampler")
## Compare chains of Gibbs sampler with MH
grid.arrange(chainsPlot_mh, chainsPlot_gibb, ncol = 2,
             top = textGrob("TracePlots of the generated Chains for the stationary Model",
                            gp = gpar(col ="#33666C",
                                      fontsize = 28, font = 2)))


param_gibbs <- apply(gibbs_statio$out.chain[,1:3], 2, mean)
param_gibbs["logsig"] <- exp(param_gibbs["logsig"] )

frame <- data.frame(Bayesian = param_gibbs, 'Frequentist(mle)' = gev_fit$mle)
row.names(frame) = c("$\\mu \\ $", "$\\sigma \\quad$", "$\\xi \\quad$")
knitr::kable(frame, align = 'l')

# Effective sample sizes :
effectiveSize(mcmc(gibbs_statio$out.chain[,1:3]))



# Function to create mcmc.lists, useful for diagnostics on chains.
'mc.listDiag3' <- function (list, subset = c("mu", "logsig", "xi")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset])
  )
}

## Predictive accuracy criterion
ic_vals <- gibbs_statio$dic.vals

# DIC Values
mc.listDiag3(gibbs_statio$out.ind)[[1]] %>%
  PissoortThesis::dic_3p(vals = ic_vals[[1]])

mc.listDiag3(gibbs_statio$out.ind)[[2]] %>%
  PissoortThesis::dic_3p(vals = ic_vals[[2]])

mc.listDiag3(gibbs_statio$out.ind)[[3]] %>%
  PissoortThesis::dic_3p(vals = ic_vals[[3]])

# WAIC Values
PissoortThesis::waic( ic_vals[[1]] )
PissoortThesis::waic( ic_vals[[2]] )
PissoortThesis::waic( ic_vals[[3]] )



#####################################################################
###########  Gibbs Sampler dealing with Nonstationarity  ############
#####################################################################


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
set.seed(101)
gibbs.trend <- PissoortThesis::gibbs.trend.own(start,
                                               propsd = c(.5, 1.9, .15, .12),
                                               iter = 1000)
colMeans(do.call(rbind, gibbs.trend$mean_acc.rates))


param.chain <- gibbs.trend$out.chain[ ,1:4]

# Effective sample sizes :
effectiveSize(mcmc(param.chain))


### Plot of the chains
PissoortThesis::chains.plotOwn(gibbs.trend$out.chain )
ggplot(gibbs.trend$out.chain) + geom_line(aes(x = iter, y = mu1)) +
  theme_piss(16,14) + labs(ylab = "mu1", xlab = "iter" )


## TracePlots
chain.mix <- cbind.data.frame(gibbs.trend$out.chain,
                              iter.chain = rep(1:500, 4))
mixchains.Own(chain.mix)
PissoortThesis::mixchains.Own(chain.mix)




library("bayesplot")
array.post <- array(unlist(gibbs.trend$out.ind), dim = c(4000/4+1, 4, 4 ),
                    dimnames = list(iterations = NULL,
                                    parameters = c("mu", "mu1", "logsig", "xi"),
                                    chains = c("chain:1", "chain:2", "chain:3",
                                               "chain:4") ))

color_scheme_set("mix-blue-red")
mcmc_trace(array.post, pars = c("mu", "mu1"),
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_trace_highlight(array.post, pars = "mu", highlight = 3)



## Gelman Coda Diagnostics
gelman.diag(mc.listDiag(gibbs.trend$out.ind), autoburnin=F)
gelman.plot(mc.listDiag(gibbs.trend$out.ind), autoburnin=F)



# Markov Chain Correlations
autocorr(mcmc(param.chain ))
autocorr.diag(mcmc(param.chain ))
autocorr.plot(mcmc(param.chain ))

crosscorr(mcmc(param.chain ))
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






## Handling of the slope parameter mu_trend
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


## Posterior Predictive Distribution
repl <- PissoortThesis::pred_post_samples()

post.pred <- apply(repl, 2, function(x) quantile(x, probs = c(0.05,0.5,0.95)))


df.postpred <- data.frame(data = max_years$data, q05 = post.pred["5%",],
                          q50 = post.pred["50%",], q95 = post.pred["95%",],
                          year = seq(1901:2016))
ggplot(df.postpred) + geom_point(aes(x = year, y = data)) +
  geom_line(aes(x = year, y = q05)) + geom_line(aes(x = year, y = q50)) +
  geom_line(aes(x = year, y = q95)) +
  ggtitle("Original data with PPD quantiles 5, 50 and 95%") + theme_piss()


## And for predictions ? 10 years here
repl2 <- PissoortThesis::pred_post_samples(n_future = 10)

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

gg1 <- PissoortThesis::Pred_Dens_ggPlot(1901, repl2)
gg2 <- PissoortThesis::Pred_Dens_ggPlot(1950, repl2)
gg3 <- PissoortThesis::Pred_Dens_ggPlot(2016, repl2)
gg4 <- PissoortThesis::Pred_Dens_ggPlot(2026, repl2)
grid.arrange(gg1, gg2, gg3, gg4, nrow = 1)


## Provide better plot with geom_joy(). (To include in package !!!!)
## Definition of parameters is straightfoward : it defines time at which we predict
# An by = is the number of densities we want to draw !
'posterior_pred_ggplot' <- function(from = 1, until = nrow(max_years$df),
                                    n_future = 0, by = 10, x_coord = c(27,35)) {

  repl2 <- PissoortThesis::pred_post_samples(from = from, until = until,
                                             n_future = n_future)

  repl2_df <- data.frame(repl2)
  colnames(repl2_df) <- seq(from + 1900, length = ncol(repl2))
  ## Compute some quantiles to later draw on the plot
  quantiles_repl2 <- apply(repl2_df, 2,
                           function(x) quantile(x , probs = c(0.05, 0.5, 0.95)) )
  quantiles_repl2 <- as.data.frame(t(quantiles_repl2))
  quantiles_repl2$year <- colnames(repl2_df)

  repl2_df_gg <- repl2_df[, seq(1, (until + n_future) - from, by = by)] %>%
    gather(year, value)

  col.quantiles <- c("5%-95%" = "cyan", "50%" = "red")
  titl <- "Posterior (predictive) densities evolution of the TX  in [1901-2021] "
  subtitl <- "Last density in 2021 is extrapolation "
  cap <- "We clearly see the linear trend from the temporal evolution of the quantiles..."
  g <- ggplot(repl2_df_gg, aes(x = value, y = as.numeric(year) )) +   # %>%rev() inside aes()
    geom_joy(aes(fill = year)) +
    geom_point(aes(x = `5%`, y = as.numeric(year), col = "5%-95%"),
               data = quantiles_repl2, size = 0.9) +
    geom_point(aes(x = `50%`, y = as.numeric(year), col = "50%"),
               data = quantiles_repl2, size = 0.9) +
    geom_point(aes(x = `95%`, y = as.numeric(year), col = "5%-95%"),
               data = quantiles_repl2, size = 0.9) +
    geom_hline(yintercept = 2016, linetype = "dashed", size = 0.3, col  = 2) +
    scale_fill_viridis(discrete = T, option = "D", direction = -1, begin = .1, end = .9) +
    scale_y_continuous(breaks = c(  seq(1901, 2016, by = by),
      seq(2016, colnames(repl2_df)[ncol(repl2_df)], by = by) ) )  +
    coord_cartesian(xlim = x_coord) +
    theme_piss(theme = theme_minimal()) +
    labs(x = "TX", y = "Year", title = titl, subtitle = subtitl, caption = cap) +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(legend.position = c(.947, .35),
          plot.subtitle = element_text(hjust = 0.5,face = "italic"),
          plot.caption = element_text(hjust = 0.1,face = "italic"))
  g
}

# Plot for the already observed data :
posterior_pred_ggplot(by = 8)

# Plot for future data (extrapolation) :
posterior_pred_ggplot(from = nrow(max_years$df), x_coord = c(30, 38),
                      n_future = nrow(max_years$df), by = 8)

## All
posterior_pred_ggplot(from = 1, x_coord = c(27, 38),
                      n_future = nrow(max_years$df), by = 12)


'mc.listDiag4' <- function (list, subset = c("mu", "mu1", "logsig", "xi")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset]),
            mcmc(list[[4]][, subset])
  )
}

## Predictive accuracy criterion
ic_vals.trend <- gibbs.trend$dic.vals

# DIC Values
mc.listDiag4(gibbs.trend$out.ind)[[1]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[1]])

mc.listDiag4(gibbs.trend$out.ind)[[2]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[2]])

mc.listDiag4(gibbs.trend$out.ind)[[3]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[3]])

mc.listDiag4(gibbs.trend$out.ind)[[4]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[4]])

# WAIC Values
PissoortThesis::waic( ic_vals.trend[[1]] )
PissoortThesis::waic( ic_vals.trend[[2]] )
PissoortThesis::waic( ic_vals.trend[[3]] )
PissoortThesis::waic( ic_vals.trend[[4]] )


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

"rl_bayes_trend" <- function(data_year, params, t = 10, m = 10 ){
    y_m <- -(1 / log(1 - 1/m))
    t <- seq(max(data_year), max(data_year) + t, 1)
    rl_m <- (params[1] + params[2] * (t-max(max_years$df$Year))) +
       (params[3] / params[4]) *  (y_m^params[4] - 1)
    g <- ggplot(data.frame(r.lvels = rl_m, years = t)) +
       geom_point(aes(x = years, y = r.lvels))
     g
  return(rl_m)
}
par_gibbs_trend[3] <- exp(par_gibbs_trend["logsig"])
rl_bayes_trend(max_years$data,
               unname(par_gibbs_trend[c("mu", "mu1", "sigma", "xi")]))

evdbayes::rl.pst(lh = "gev")





## Make the time component in squared ? ie mu = mu_0 + mu_1 * t + mu_2 * t^2
##############################################################################


fn2 <- function(par, data) -log_post2(par[1], par[2], par[3],
                                      par[4], par[5], data )
param2 <- c( mu0 = mean(max_years$df$Max), mu1 = 0, mu2 = 0,
             logsig = log(sd(max_years$df$Max)), xi =  -0.1 )
opt2 <- optim(param2, fn2, data = max_years$data,
              method = "BFGS", hessian = T)
opt2

# Starting Values
set.seed(100)
start <- list() ;   k <- 1
while(k < 6) { # starting value is randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt2$par, solve(opt2$hessian)))
  svlp <- log_post2(sv[1], sv[2], sv[3], sv[4], sv[5], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start[[k]] <- sv
    k <- k + 1
  }
}
# k chains with k different starting values
set.seed(101)
gibbs.trend2 <- gibbs.trend2.own(start, propsd = c(.5, 1.4, 3.5, .2, .15),
                                 iter = 1000)
colMeans(do.call(rbind, gibbs.trend2$mean_acc.rates))

param.chain2 <- gibbs.trend2$out.chain[ ,1:5]


### Plot of the chains
chains.plotOwn(gibbs.trend2$out.chain)
grid.arrange(
  ggplot(gibbs.trend2$out.chain) + geom_line(aes(x = iter, y = mu1)) +
    theme_piss(16,14) + labs(ylab = "mu1", xlab = "iter"),
  ggplot(gibbs.trend2$out.chain) + geom_line(aes(x = iter, y = mu2)) +
    theme_piss(16,14) + labs(ylab = "mu2", xlab = "iter")    )


## TracePlots
chain.mix <- cbind.data.frame(gibbs.trend2$out.chain,
                              iter.chain = rep(1:500, 4))
mixchains.Own(chain.mix)
ggplot(chain.mix, aes(x = iter.chain, y = mu1, col = as.factor(chain.nbr))) +
  geom_line() + theme_piss(18,16, theme_classic()) +
  scale_colour_brewer(name = "chain nr", palette = "Set1") +
  guides(colour = guide_legend(override.aes = list(size= 1.2)))
ggplot(chain.mix, aes(x = iter.chain, y = mu2, col = as.factor(chain.nbr))) +
  geom_line() + theme_piss(18,16, theme_classic()) +
  scale_colour_brewer(name = "chain nr", palette = "Set1") +
  guides(colour = guide_legend(override.aes = list(size= 1.2)))



'mc.listDiag5' <- function (list, subset = c("mu", "mu1", "mu2",
                                             "logsig", "xi")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset]),
            mcmc(list[[4]][, subset]),
            mcmc(list[[5]][, subset])
  )
}
## Predictive accuracy criterion

ic_vals.trend2 <- gibbs.trend2$dic.vals

# DIC Values. "1:5" takes the 5 parameters of the model, see function
mc.listDiag5(gibbs.trend2$out.ind)[[1]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[1]])

mc.listDiag5(gibbs.trend2$out.ind)[[2]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[2]])

mc.listDiag5(gibbs.trend2$out.ind)[[3]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[3]])

mc.listDiag5(gibbs.trend2$out.ind)[[4]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[4]])

mc.listDiag5(gibbs.trend2$out.ind)[[5]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[5]])


# WAIC Values
PissoortThesis::waic( ic_vals.trend2[[1]] )
PissoortThesis::waic( ic_vals.trend2[[2]] )
PissoortThesis::waic( ic_vals.trend2[[3]] )
PissoortThesis::waic( ic_vals.trend2[[4]] )
PissoortThesis::waic( ic_vals.trend2[[5]] )
## Interestingly here, all the criterion are lower than for simple trend
# and this models should be preferred... different result from frequentist !!


## Cross-validation




######## Model with inear trend + varying scale parameter (exp link)  ####
########   mu = mu0 + mu1 * t,     logsig = sig0 + sig1 * tt     #########
##########################################################################


fn3 <- function(par, data) -log_post3(par[1], par[2], par[3],
                                      par[4], par[5], data )
param3 <- c( mu0 = mean(max_years$df$Max), mu1 = 0,
             sig0 = log(sd(max_years$df$Max)), sig1 = 0,  xi =  -0.1 )
opt3 <- optim(param3, fn3, data = max_years$data,
              method = "BFGS", hessian = T)
opt3

# Starting Values
set.seed(100)
start <- list() ; k <- 1
while(k < 6) { # starting value is randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt3$par, solve(opt3$hessian)))
  svlp <- log_post1(sv[1], sv[2], sv[3], sv[4], sv[5], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start[[k]] <- sv
    k <- k + 1
  }
}
# k chains with k different starting values
set.seed(101)
gibbs.trend.sig3 <- gibbs.trend.sig3own(start,
                                        propsd = c(.45, 1.3, .2, .5, .1),
                                        iter = 1e3)
colMeans(do.call(rbind, gibbs.trend.sig3$mean_acc.rates))

param.chain3 <- gibbs.trend.sig3$out.chain[, 1:5]


### Plot of the chains
grid.arrange(
  ggplot(gibbs.trend.sig3$out.chain) + geom_line(aes(x = iter, y = mu)) +
    theme_piss(16,14) + labs(ylab = "mu", xlab = "iter"),
  ggplot(gibbs.trend.sig3$out.chain) + geom_line(aes(x = iter, y = mu1)) +
    theme_piss(16,14) + labs(ylab = "mu1", xlab = "iter"),
  ggplot(gibbs.trend.sig3$out.chain) + geom_line(aes(x = iter, y = sig0)) +
    theme_piss(16,14) + labs(ylab = "sig0", xlab = "iter"),
  ggplot(gibbs.trend.sig3$out.chain) + geom_line(aes(x = iter, y = sig1)) +
    theme_piss(16,14) + labs(ylab = "sig1", xlab = "iter"),
  ggplot(gibbs.trend.sig3$out.chain) + geom_line(aes(x = iter, y = xi)) +
    theme_piss(16,14) + labs(ylab = "xi", xlab = "iter")    )


## TracePlots
chain.mix <- cbind.data.frame(gibbs.trend2$out.chain,
                              iter.chain = rep(1:500, 4))
mixchains.Own(chain.mix)
ggplot(chain.mix, aes(x = iter.chain, y = mu1, col = as.factor(chain.nbr))) +
  geom_line() + theme_piss(18,16, theme_classic()) +
  scale_colour_brewer(name = "chain nr", palette = "Set1") +
  guides(colour = guide_legend(override.aes = list(size= 1.2)))




## Predictive accuracy criterion

ic_vals.trend.sc <- gibbs.trend.sig3$dic.vals

# DIC Values. "1:5" takes the 5 parameters of the model, see function

mc.listDiag5(gibbs.trend.sig3$out.ind, 1:5)[[1]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[1]])

mc.listDiag5(gibbs.trend.sig3$out.ind, 1:5)[[2]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[2]])

mc.listDiag5(gibbs.trend.sig3$out.ind, 1:5)[[3]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[3]])

mc.listDiag5(gibbs.trend.sig3$out.ind, 1:5)[[4]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[4]])

mc.listDiag5(gibbs.trend.sig3$out.ind, 1:5)[[5]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[5]])

# WAIC Values
waic( ic_vals.trend.sc[[1]] )
waic( ic_vals.trend.sc[[2]] )
waic( ic_vals.trend.sc[[3]] )
waic( ic_vals.trend.sc[[4]] )
waic( ic_vals.trend.sc[[5]] )


# Comparing These values with the ones obtained with simple linear trend, all are
# again for the complex model .....



##################### Comparisons  ##########################################
#############################################################################
library(bbefkr)  # Compute the marginal likelihood using Chib's (1995) method
# To be able to compute the Bayes factor.
logdensity_admkr()




## Cross-validation

library(loo)

extract_log_lik()



#save.image("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1_bayes.Rdata")

#### See end of bayesian01 for predictive dist, quantiles, and other...
