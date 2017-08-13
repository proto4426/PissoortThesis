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
library(reshape2)

library(PissoortThesis)



## As we aim at building models sequentially,  we will begin by the Gumbel model
## It is models the lowest number of degrees of freedom.


## MEtropolis-Hastlings for the stationary Gumbel model ###
###########################################################


# Optimize Posterior Density Function to find starting values
fn <- function(par, data) -log_post_gumb(par[1], par[2], data)
param <- c(mean(max_years$df$Max), log(sd(max_years$df$Max)))
opt0 <- nlm(fn, param, data = max_years$data,
           hessian=T, iterlim = 1e5)
start0 <- opt0$estimate  ;     names(start0) <- c("mu", "logsig")
Sig <- solve(opt0$hessian)
ev <- eigen( (2.4/sqrt(2))^2 * Sig)
varmat0 <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)

set.seed(101)
iter <- 2e3
mh.mcmc_gumb <- PissoortThesis::MH_mcmc.own(start0,
                                            varmat0 %*% c(.9, 1.1),
                                            iter, burnin = iter/4)
# Display the mean acceptance rate
mh.mcmc_gumb$mean.acc_rates



##### Gumbel Model with the Gibbs Sampler
#########################################

set.seed(101)
iter <- 2e3
gib.mcmc_gumbel <- PissoortThesis::gibbs_mcmc.own(start0, iter = iter,
                                                  nbr.chain = 1,
                                                  Gumbel = T,
                                                  propsd = c(.42, .12),
                                                  burnin = iter/4)


#########   MEtropolis-Hastlings for the stationary GEV model #######
#####################################################################


# Optimize Posterior Density Function to find starting values
fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
param <- c(mean(max_years$df$Max), log(sd(max_years$df$Max)), 0.1 )
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
# Display the mean acceptance rate
mh.mcmc1$mean.acc_rates

param_mean_mh <- apply(mh.mcmc1$out.chain[,1:3], 2, mean)
chainsPlot_mh <- chains.plotOwn(mh.mcmc1$out.chain,
                                post.mean.green = param_mean_mh,
                                title = "Using Metropolis-Hastings Algorithm")


## Burn in 1/4 of values  :
mh.mcmc.out <- mh.mcmc1$out.chain[-(1:iter/4), ]

# Effective sample sizes :
effectiveSize(mcmc(mh.mcmc.out[,1:3]))




##########  GIBBS sampler for the stationary GEV model  #######
##############################################################


## Advised value for the proposal variance
#prop_var <- sqrt( (2.4/sqrt(1))^2 * solve(opt$hessian) )


set.seed(100)
iter <- 2e3
# Same starting point as for MH. We also keep only one chain
# i.e. no different starting values so far.
gibbs_statio <- PissoortThesis::gibbs_mcmc.own(start, iter = iter,
                                               nbr.chain = 1,
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

## Compare the parameters from Gibbs and MH
param_gibbs <- apply(gibbs_statio$out.chain[,1:3], 2, mean)

# Transform back the scale parameter
param_mean_mh["sigma"] <- exp(param_mean_mh["logsig"] )
#param_mean_mh <- param_mean_mh[names(param_mean_mh) %in% "logsig" == F]
param_gibbs["sigma"] <- exp(param_gibbs["logsig"] )
#param_gibbs <- param_gibbs[names(param_gibbs) %in% "logsig" == F]

frame <- data.frame(Bayesian.mh = param_mean_mh[c("mu", "sigma", "xi")],
                    Bayesian.gib = param_gibbs[c("mu", "sigma", "xi")],
                    'Frequentist(mle)' = gev_fit$mle)
row.names(frame) = c("$\\mu \\ $", "$\\sigma \\quad$", "$\\xi \\quad$")
#knitr::kable(frame, align = 'l')
stargazer::stargazer(frame, summary = F)


# Effective sample sizes :
effectiveSize(mcmc(gibbs_statio$out.chain[,1:3]))



# Function to create mcmc.lists, useful for diagnostics on chains.
'mc.listDiag3' <- function (list, subset = c("mu0", "logsig", "xi")) {
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

# WAIC Values
PissoortThesis::waic( ic_vals[[1]] )




#####################################################################
###########  Gibbs Sampler dealing with Nonstationarity : ###########
##### Linear model on the location :  mu(t) = mu0 + mu1 * t  ########
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
set.seed(104)
start <- list() ;   k <- 1
while(k < 5) { # starting values are randomly selected from a distribution
               # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt$par, 50 * solve(opt$hessian)))
  svlp <- log_post1(sv[1], sv[2], sv[3], sv[4], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start[[k]] <- sv
    k <- k + 1
  }
}
mat_startvalues <- matrix(unlist(start), ncol = 4, byrow = T,
                          dimnames = list(c("start[[1]]", "start[[2]]",
                                            "start[[3]]", "start[[4]]"),
                                          c("$\\mu_0 \\ $","$\\mu_{trend}$",
                                            "$\\logsig$", "$\\xi$")))
knitr::kable(mat_startvalues, align = "c")
stargazer::stargazer(mat_startvalues, summary = F)

# k chains with k different starting values. Number of iter is for each chain.
set.seed(102)
iter.by.chain <- 1000
gibbs.trend <- PissoortThesis::gibbs.trend.own(start,
                                               propsd = c(.5, 1.9, .15, .12),
                                               iter = iter.by.chain,
                                               burnin = ceiling(iter.by.chain/2 + 1))
gibbs.trend$mean_acc.rates
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
                              iter.chain = rep(501:1000, 4))
chain_mix_gg <- mixchains.Own(chain.mix, burnin = 501)
title = "TracePlots of the generated Chains "
PissoortThesis::grid_arrange_legend(chain_mix_gg$gmu, chain_mix_gg$gmutrend,
                                     ncol = 2,
                                    top = grid::textGrob(title,
                                      gp = grid::gpar(col = "black",
                                                      fontsize = 25, font = 4)) )
grid.arrange(chain_mix_gg$glogsig, chain_mix_gg$gxi, ncol = 2)

library("bayesplot")
array.post <- array(unlist(gibbs.trend$out.ind), dim = c(4000/4+1, 4, 4 ),
                    dimnames = list(iterations = NULL,
                                    parameters = c("mu0", "mu1", "logsig", "xi"),
                                    chains = c("chain:1", "chain:2", "chain:3",
                                               "chain:4") ))
color_scheme_set("mix-blue-red")
mcmc_trace(array.post,
           facet_args = list(ncol = 1, strip.position = "left"))
mcmc_trace_highlight(array.post, highlight = 3)



'mc.listDiag4' <- function (list, subset = c("mu0", "mu1", "logsig", "xi")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset]),
            mcmc(list[[4]][, subset])
  )
}

## Gelman Coda Diagnostics
gelman.diag(mc.listDiag4(gibbs.trend$out.ind), autoburnin=F)
gelman.plot(mc.listDiag4(gibbs.trend$out.ind), autoburnin=F, auto.layout = T)
# In ggplot : Put all on the same y-scales
gp.dat <- gelman.plot(mc.listDiag4(gibbs.trend$out.ind), autoburnin=F)
df = data.frame(bind_rows(as.data.frame(gp.dat[["shrink"]][,,1]),
                          as.data.frame(gp.dat[["shrink"]][,,2])),
                q = rep(dimnames(gp.dat[["shrink"]])[[3]],
                      each = nrow(gp.dat[["shrink"]][,,1])),
                last.iter = rep(gp.dat[["last.iter"]], length(gp.dat)))
df_gg <-melt(df, c("q","last.iter"), value.name = "shrink_factor")
ggplot(df_gg,
       aes(last.iter, shrink_factor, colour=q, linetype=q)) +
  geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
  geom_line() +
  geom_hline(yintercept = 1.1, colour = "green4", linetype = "dashed", size = 0.3) +
  scale_y_continuous(breaks = c(1, 1.1, 1.5, 2, 3, 4 ),
                     labels = c(1, 1.1, 1.5, 2, 3, 4 )) +
  #ggtitle("Gelman Rubin dignostic : R-hat Statistic") +
  facet_wrap(~variable,
             labeller= labeller(.cols=function(x) gsub("V", "Chain ", x))) +
  labs(x="Last Iteration in Chain", y="Shrink Factor",
       colour="Quantile", linetype="Quantile",
       subtitle = "Gelman Rubin diagnostic : R-hat Statistic") +
  scale_linetype_manual(values=c(2,1)) +
  theme_piss() +
  theme(strip.text = element_text(size=15),
        plot.subtitle = element_text(size = 21, hjust = 0.5,
                                     colour = "#33666C", face = "bold"))


# Transform back the scale parameter
param.chain[, "sigma"] <- exp(param.chain[,"logsig"])

# Markov Chain Correlations
autocorr(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")] ))
autocorr.diag(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")] ))
autocorr.plot(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")]  ))

crosscorr.plot(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")]  ),
                title = "Cross-correlation")
title("Cross-correlation")
# In ggplot
library(ggcorrplot)
ggcorrplot(crosscorr(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")])),
           hc.order = TRUE, type = "lower", lab = TRUE, title = "Cross-correlation",
           ggtheme = PissoortThesis::theme_piss)
# Compare it with Fisher information matrix ('frequentist') by MLE with ismev
cr.corr_lin <- crosscorr(mcmc(param.chain[, c("mu0", "mu1", "sigma", "xi")]))
dimnames(gev_nonstatio$cov) <- dimnames(cr.corr_lin)
# Transform it to correlation for comparison
cov2cor(gev_nonstatio$cov)  ;    cr.corr_lin

## Geweke
geweke <- geweke.diag(mcmc(param.chain))
2*dnorm(geweke$z)
geweke.plot(mcmc(param.chain), nbins = 20)


## Raftery Coda Diagnostics
# For each chain individually
raf1 <- raftery.diag(mc.listDiag4(gibbs.trend$out.ind)[[1]], q=0.05, r=0.02, s=0.95)
raf2 <- raftery.diag(mc.listDiag4(gibbs.trend$out.ind)[[2]], q=0.05, r=0.02, s=0.95)
raf3 <- raftery.diag(mc.listDiag4(gibbs.trend$out.ind)[[3]], q=0.05, r=0.02, s=0.95)
raf4 <- raftery.diag(mc.listDiag4(gibbs.trend$out.ind)[[4]], q=0.05, r=0.02, s=0.95)
set.seed(12)
raf <- sample(list(raf1$resmatrix, raf2$resmatrix, raf3$resmatrix, raf4$resmatrix), 1)
stargazer::stargazer(raf, summary = F)
# For the complete chains
raf_tot <- raftery.diag(mcmc(gibbs.trend$out.chain[, c("mu0", "mu1", "logsig", "xi")]),
                     q=0.05, r=0.02, s=0.95)
stargazer::stargazer(raf_tot$resmatrix, summary = F)



## Recompute the chain by increasing the number of iterations
set.seed(102)
iter.by.chain2 <- 2500
gibbs.trend_more <- PissoortThesis::gibbs.trend.own(start,
                                               propsd = c(.5, 1.9, .15, .12),
                                               iter = iter.by.chain2,
                                               burnin = ceiling(iter.by.chain2/4 + 1))
gibbs.trend_more$mean_acc.rates
colMeans(do.call(rbind, gibbs.trend_more$mean_acc.rates))


param.chain_more <- gibbs.trend_more$out.chain[ ,1:4]

# Effective sample sizes :
effectiveSize(mcmc(param.chain_more))




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
par_gibbs_trend <- apply(gibbs.trend$out.chain[,c("mu0", "mu1", "logsig", "xi")],
                         2, mean) # Average over the (3) generated chains
par_gibbs_trend["sigma"] <- exp(par_gibbs_trend["logsig"] )
names(par_gibbs_trend["logsig"]) <- "sigma"
frame <- data.frame(Bayesian = par_gibbs_trend,
                    'Frequentist(mle)' = gev_nonstatio$mle)
row.names(frame) = c("$\\mu_0 \\ $","$\\boxed{\\mu_1} \\quad$",
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
g1 <- mcmc_dens(param.chain, pars = c("mu0")) +
  geom_vline(xintercept = tab_quantiles['mu0', "2.5%"], col = "red") +
  geom_vline(xintercept = tab_quantiles['mu0', "97.5%"], col = "red")
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

mut <- par_gibbs_trend["mu0"] + par_gibbs_trend["mu1"] * t_bayes
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


## Provide better plot with geom_joy(). (can be included in the package !)
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
  cap <- "We clearly see the linear trend with the temporal evolution of the quantiles...
          And the uncertainty from their deviation  from the median"

  g <- ggplot(repl2_df_gg, aes(x = value, y = as.numeric(year) )) +  # %>%rev() inside aes()
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
rl.post_gg(gibb1$out.chain[, c("mu0", "sigma", "xi")], method = "gev")
rl.pred_gg(gibb1$out.chain[, c("mu0", "sigma", "xi")], method = "gev",
           qlim = c(30, 40), period = c(1, 5, 15))



# Nonstationary (linear trend)
library(evdbayes)
gibbs.trend$out.chain[,"sigma"] <- exp(gibbs.trend$out.chain[,"logsig"])
rl.data <-gibbs.trend$out.chain[, c("mu0", "mu1", "sigma", "xi")]
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
               unname(par_gibbs_trend[c("mu0", "mu1", "sigma", "xi")]))

evdbayes::rl.pst(lh = "gev")





## Make the time component in squared ? ie mu(t) = mu_0 + mu_1 * t + mu_2 * t^2
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



'mc.listDiag5' <- function (list, subset = c("mu0", "mu1", "mu2",
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
waic.trend2.1 <- PissoortThesis::waic( ic_vals.trend2[[1]] )
waic.trend2.2 <- PissoortThesis::waic( ic_vals.trend2[[2]] )
waic.trend2.3 <- PissoortThesis::waic( ic_vals.trend2[[3]] )
waic.trend2.4 <- PissoortThesis::waic( ic_vals.trend2[[4]] )
waic.trend2.5 <- PissoortThesis::waic( ic_vals.trend2[[5]] )
cat(waic.trend2.1, waic.trend2.2, waic.trend2.3, waic.trend2.4, waic.trend2.5)
## Interestingly here, all the criterion are lower than for simple trend
# and this models should be preferred... different result from frequentist !!

waic.mean.trend2 <- mean(waic.trend2.1, waic.trend2.2, waic.trend2.3,
                       waic.trend2.4, waic.trend2.5)

## Cross-validation




######## Model with inear trend + varying scale parameter (exp link)  ####
########   mu(t) = mu0 + mu1 * t,     logsig(t) = sig0 + sig1 * tt     #########
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
    theme_piss(16,14) + labs(ylab = "mu0", xlab = "iter"),
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




#### Predictive accuracy criterion

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
waic.sig3.1 <- PissoortThesis::waic( ic_vals.trend.sc[[1]] )
waic.sig3.2 <- PissoortThesis::waic( ic_vals.trend.sc[[2]] )
waic.sig3.3 <- PissoortThesis::waic( ic_vals.trend.sc[[3]] )
waic.sig3.4 <- PissoortThesis::waic( ic_vals.trend.sc[[4]] )
waic.sig3.5 <- PissoortThesis::waic( ic_vals.trend.sc[[5]] )
cat(waic.sig3.1, waic.sig3.2, waic.sig3.3, waic.sig3.4, waic.sig3.5)

waic.mean.sig3 <- mean(waic.sig3.1, waic.sig3.2,
                       waic.sig3.3, waic.sig3.4, waic.sig3.5)
# Comparing These values with the ones obtained with simple linear trend, all are
# again for the complex model .....



##################### Comparisons  ##########################################
#############################################################################
library(loo)


comp.waic.statio <- waic(ic_vals[[1]])

comp.waic.trend <- waic(rbind(ic_vals.trend[[1]], ic_vals.trend[[2]],
                               ic_vals.trend[[3]], ic_vals.trend[[4]]))

comp.waic.trend2 <- waic(rbind(ic_vals.trend2[[1]], ic_vals.trend2[[2]],
                               ic_vals.trend2[[3]], ic_vals.trend2[[4]],
                               ic_vals.trend2[[5]]))
loo(rbind(ic_vals.trend2[[1]], ic_vals.trend2[[2]],
           ic_vals.trend2[[3]], ic_vals.trend2[[4]],
           ic_vals.trend2[[5]]))

comp.waic.sig3 <- waic(rbind(ic_vals.trend.sc[[1]], ic_vals.trend.sc[[2]],
                             ic_vals.trend.sc[[3]], ic_vals.trend.sc[[4]],
                             ic_vals.trend.sc[[5]]))

cat(comp.waic.statio$waic, comp.waic.trend$waic,
    comp.waic.trend2$waic, comp.waic.sig3$waic)
loo::compare(comp.waic.statio, comp.waic.trend,
             comp.waic.trend2, comp.waic.sig3)



library(bbefkr)  # Compute the marginal likelihood using Chib's (1995) method
# To be able to compute the Bayes factor.
logdensity_admkr()

library(BayesFactor)
data(puzzles)
bf = anovaBF(RT ~ shape*color + ID, whichRandom = "ID", data = puzzles)
bf


## Cross-validation


# Read in and prepare the data
wells <- read.csv("wells.csv")
N <- nrow(wells)
X <- cbind(rep(1,N), wells$dist100, wells$arsenic)
y <- wells$y
P <- ncol(X)
# Fit the model with Stan
fit_1 <- stan("logistic.stan")
print(fit_1, "b")
# Compute LOO
log_lik_1 <- extract_log_lik(fit_1)
loo_1 <- loo(log_lik_1)
print(loo_1)


#save.image("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1_bayes.Rdata")

#### See end of bayesian01 for predictive dist, quantiles, and other...
