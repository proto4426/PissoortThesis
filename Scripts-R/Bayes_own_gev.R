#setwd('/home/proto4426/Documents/Extreme/R resources/IRM')
#load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1.Rdata")
load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1_bayes.Rdata")


# Bayesian Analysis constrcuted with functions in package PissoortThesis
######################################################################
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
data("max_years")



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


# Starting Values
set.seed(10)
start0 <- list() ;   k <- 1
while(k < 5) { # starting values are randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt0$estimate, 50 * solve(opt0$hessian)))
  svlp <- log_post_gumb(sv[1], sv[2], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start0[[k]] <- sv
    k <- k + 1
  }
}

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
iter <- 1e3
gib.mcmc_gumbel <- gibbs_mcmc.own(start0, iter = iter,
                                                  nbr.chain = 4,
                                                  Gumbel = T,
                                                  propsd = c(.42, .12),
                                                  burnin = iter/4)
## (For  computation time comparisons we take these values for the parameter here)




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

# Starting Values
set.seed(101)
start <- list() ;   k <- 1
while(k < 5) { # starting values are randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt$estimate, 50 * solve(opt$hessian)))
  svlp <- log_post0(sv[1], sv[2], sv[3], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start[[k]] <- sv
    k <- k + 1
  }
}


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
iter <- 1e3
# Same starting point as for MH. We also keep only one chain
# i.e. no different starting values so far.
gibbs_statio <- gibbs_mcmc.own(start, iter = iter,
                                               nbr.chain = 4,
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




#####################################################################
###########  Gibbs Sampler dealing with Nonstationarity : ###########
##### Linear model on the location :  mu(t) = mu0 + mu1 * t  ########
#####################################################################

data <- max_years$data


fn <- function(par, data) -log_post1(par[1], par[2], par[3],
                                     par[4],rescale.time = T, data)
param <- c(mean(max_years$df$Max), 0, log(sd(max_years$df$Max)), -0.1 )
opt <- optim(param, fn, data = max_years$data,
             method = "BFGS", hessian = T)
# opt <- nlm(fn, param, data = max_years$data,
#     hessian=T, iterlim = 1e5)
opt

# Starting Values
set.seed(10)
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
                                               burnin = ceiling(iter.by.chain/2 + 1),
                                               keep.same.seed = 12)
gibbs.trend$mean_acc.rates
colMeans(do.call(rbind, gibbs.trend$mean_acc.rates))


head( param.chain <- gibbs.trend$out.chain[ ,1:4] )

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
                                      gp = grid::gpar(col = "#33666C",
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


# Function to create mcmc.lists, useful for diagnostics on chains.
'mc.listDiag4' <- function (list, subset = c("mu0", "mu1", "logsig", "xi")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset]),
            mcmc(list[[4]][, subset])
  )
}

## Gelman Coda Diagnostics
Rhat <- gelman.diag(mc.listDiag4(gibbs.trend$out.ind), autoburnin=F)
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



## Recompute the chain by decreasing the Burn-in period
set.seed(102)
# keep.same.seed  sets a seed at each iterations that is the integer you specify
# times the iteration number.
gibbs.trend <- PissoortThesis::gibbs.trend.own(start,
                                               propsd = c(.5, 0.008, .15, .12), # 0.008 ; 1.9
                                               iter = iter.by.chain,
                                               burnin = ceiling(iter.by.chain2/4 + 1),
                                               keep.same.seed = 123,
                                               rescale.time = F)
gibbs.trend$mean_acc.rates
colMeans(do.call(rbind, gibbs.trend$mean_acc.rates))
head(param.chain <- gibbs.trend$out.chain[ ,1:4])
# Effective sample sizes :
effectiveSize(mcmc(param.chain))



##### Inference  ########

par_gibbs_trend <- apply(gibbs.trend$out.chain[,c("mu0", "mu1", "logsig", "xi")],
                         2, mean) # Average over the (4) generated chains

## Handling of the slope parameter mu_trend
t_bayes <- round(( min(max_years$df$Year):max(max_years$df$Year) -
                     mean(max_years$df$Year) ) / length(max_years$data),4)
t_freq <- seq(1, length(max_years$data),1)

mut <- par_gibbs_trend["mu0"] + par_gibbs_trend["mu1"] * t_bayes
mu_freq <- gev_nonstatio$mle[1] + gev_nonstatio$mle[2] * t_freq
# See the values of mu(t) in Bayesian VS in frequentist
pander(t(matrix(c(mut[1:4], mut[(length(mut)-3):length(mut)],
                  mu_freq[1:4], mu_freq[(length(mu_freq)-3):length(mu_freq)]),
                ncol = 2, dimnames = list(c(rep("First values", 4),
                                            rep("Last values", 4)),
                                          c("Frequentist",  "Bayesian")))))
( mu_freq[length(mu_freq)] - mu_freq[1] ) / length(mu_freq)
( mut[length(mut)] - mut[1] ) / length(mut)


## Summary And Parameter Table
param.chain$sigma <- exp(param.chain$logsig)
tab_quantiles <- as.data.frame(summary(mcmc(param.chain))$quantiles)

stargazer::stargazer(tab_quantiles, summary = F)

mean.post <- apply(param.chain, 2, mean)

# Transform back mu_1
mut <- tab_quantiles$`2.5%`[1] + tab_quantiles$`2.5%`[2] * t_bayes
cat("q.2.5% is ", q2.5_mu1Trans <- ( mut[length(mut)] - mut[1] )  / length(mut) )
mut <- tab_quantiles$`25%`[1] + tab_quantiles$`25%`[2] * t_bayes
cat("q.25% is ", q25_mu1Trans <- ( mut[length(mut)] - mut[1] )  / length(mut) )
mut <- tab_quantiles$`50%`[1] + tab_quantiles$`50%`[2] * t_bayes
cat("q.50% is ", q50_mu1Trans <- ( mut[length(mut)] - mut[1] )  / length(mut) )
mut <- tab_quantiles$`75%`[1] + tab_quantiles$`75%`[2] * t_bayes
cat("q.75% is ", q75_mu1Trans <- ( mut[length(mut)] - mut[1] )  / length(mut) )
mut <- tab_quantiles$`97.5%`[1] + tab_quantiles$`97.5%`[2] * t_bayes
cat("q.97.5% is ", q975_mu1Trans <- ( mut[length(mut)] - mut[1] )  / length(mut) )


## Comparisons with Frequentist's results (GEV)
par_gibbs_trend["sigma"] <- exp(par_gibbs_trend["logsig"] )
names(par_gibbs_trend["logsig"]) <- "sigma"
frame <- data.frame(Bayesian = par_gibbs_trend[c("mu0", "mu1", "sigma", "xi")],
                    'Frequentist(mle)' = gev_nonstatio$mle)
knitr::kable(frame, align = 'l')


## Credible intervals (quantile-based) see above quantiles ! Densities are quite
# symetric, so it is not a bad idea to use this method.

## HPD intervals
library(HDInterval)
hpd_mu0 <- hdi(param.chain$mu0)
hpd_mu1 <- hdi(param.chain$mu1)
hpd_logsigma <- hdi(param.chain$logsig)
hpd_xi <- hdi(param.chain$xi)
hpd_sigma <- hdi(param.chain$sigma)
hpd95 <- data.frame(mu0 = c(hpd_mu0), mu1 = c(hpd_mu1),
                    logsig = c(hpd_logsigma),
                    xi = c(hpd_xi), sig = c(hpd_sigma))

hpd_mu0.75 <- hdi(param.chain$mu0,credMass = 0.5)
hpd_logsigma.75 <- hdi(param.chain$logsig, credMass = 0.5)
hpd_sigma.75 <- hdi(param.chain$sigma, credMass = 0.5)
hpd_xi.75 <- hdi(param.chain$xi, credMass = 0.5)
hpd_mu1.75 <- hdi(param.chain$mu1, credMass = 0.5)
hpd75 <- data.frame(mu0 = c(hpd_mu0.75), mu1 = c(hpd_mu1.75),
                    logsig = c(hpd_logsigma.75),
                    xi = c(hpd_xi.75), sig = c(hpd_sigma.75))
### Comparisons of all the ci :

## Densities of the parameters with their quantile-based and  HPD 0.95 intervals
color_scheme_set("brightblue")
col.intervals <- c("Quantile" = "red", "HPD" = "green")

'legend.things'  <-
  list(scale_color_manual(name = "Intervals", values = col.intervals),
    theme_piss(legend.position = c(0.92, 0.5)),
    theme(legend.background = element_rect(colour = "transparent",size = 0.5))
   )

g1 <- mcmc_dens(param.chain, pars = c("mu0")) +
  geom_vline(aes(xintercept = tab_quantiles['mu0', "2.5%"],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = tab_quantiles['mu0', "97.5%"],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_mu0[[1]],
             col = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_mu0[[2]],
             col = "HPD"), linetype = "dashed") +
  legend.things
g2 <- mcmc_dens(param.chain, pars = c("logsig")) +
  geom_vline(aes(xintercept = tab_quantiles['logsig', "2.5%"],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = tab_quantiles['logsig', "97.5%"],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_logsigma[[1]],
             col = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_logsigma[[2]],
             col = "HPD"), linetype = "dashed") +
  legend.things
g3 <- mcmc_dens(param.chain, pars = c("xi"))+
  geom_vline(aes(xintercept = tab_quantiles['xi', '2.5%'],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = tab_quantiles['xi', "97.5%"],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_xi[[1]],
             col = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_xi[[2]],
             col = "HPD"), linetype = "dashed")+
  legend.things
g4 <- mcmc_dens(param.chain, pars = c("mu1")) +
  geom_vline(aes(xintercept = tab_quantiles['mu1', '2.5%'],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = tab_quantiles['mu1', "97.5%"],
             col = "Quantile"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_mu1[[1]],
             col = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = hpd_mu1[[2]],
             col = "HPD"), linetype = "dashed") +
  legend.things

title <- "Posterior densities of the parameters and Bayesian intervals"
grid.arrange( g1, g4, g2, g3, nrow = 2,
              top = grid::textGrob(title,
                                   gp = grid::gpar(col = "#33666C",
                                                   fontsize = 25,
                                                   font = 4), vjust = 0.4))


## Compare the frequentists with the Bayesian intervals

gev_freq_mu0_prof75 <- ci(gev_nonstatio_prof, method="proflik",verbose=T,
                          type = "parameter", which.par = 1, nint = 1000, alpha = 0.25)
gev_freq_mu1_prof75 <- ci(gev_nonstatio_prof, method="proflik",verbose=T,
                          type = "parameter", which.par = 2, nint = 1000, alpha = 0.25)
gev_freq_sig_prof75 <- ci(gev_nonstatio_prof, method="proflik",verbose=T,
                          type = "parameter", which.par = 3, nint = 1000, alpha = 0.25)
gev_freq_xi_prof75 <- ci(gev_nonstatio_prof, method="proflik",verbose=T,
                         type = "parameter", which.par = 4, nint = 1000, alpha = 0.25)

# Function to compute the intervals from porfile likelihood and from bayesian quantiles
'intervals_compareByParam' <- function(proflik.ci95 = gev_freq_xi_prof95,
                                      proflik.ci75 = gev_freq_xi_prof75,
                                      which = 4) {

  int_freq.norm <- data.frame(ll = gev_freq_norm95[which, 1],
                              l = gev_freq_norm75[which, 1],
                              m = gev_nonstatio$mle[which],
                              h = gev_freq_norm75[which, 3],
                              hh = gev_freq_norm95[which, 3] )
  int_freq.prof <- data.frame(ll = proflik.ci95[[1]],
                            l = proflik.ci75[[1]],
                            m = gev_nonstatio$mle[which],
                            h = proflik.ci75[[3]],
                            hh = proflik.ci95[[3]] )
# browser()
  # hanle the transformation of sigma
  if(which == 3)   which <- which + 2   #

  ## The  If condition to handle the rescaling problem of mu_1
  if(which == 2 ) {
    # int_bayes_quant <- data.frame(ll = q2.5_mu1Trans,
    #                                        l = q25_mu1Trans,
    #                                        m = q50_mu1Trans,
    #                                        h = q75_mu1Trans,
    #                                        hh = q975_mu1Trans)
    int_bayes_quant <- data.frame(ll = tab_quantiles$`2.5%`[which],
                                  l = tab_quantiles$`25%`[which],
                                  m = tab_quantiles$`50%`[which],
                                  h = tab_quantiles$`75%`[which],
                                  hh = tab_quantiles$`97.5%`[which])
  int_bayes_hpd <- data.frame(ll = hpd_mu1.trans[[1]],
                              l = hpd_mu1_trans.75[[1]],
                              m = q50_mu1Trans,
                              h = hpd_mu1_trans.75[[2]],
                              hh = hpd_mu1.trans[[2]])
  }
  else  {  int_bayes_quant <- data.frame(ll = tab_quantiles$`2.5%`[which],
                             l = tab_quantiles$`25%`[which],
                             m = tab_quantiles$`50%`[which],
                             h = tab_quantiles$`75%`[which],
                             hh = tab_quantiles$`97.5%`[which])

  int_bayes_hpd <- data.frame(ll = hpd95[1,which],
                                l = hpd75[1, which],
                                m = tab_quantiles$`50%`[which],
                                h = hpd75[2, which],
                                hh = hpd95[2, which])
  }

  df <- cbind.data.frame("Frequentist.Normal" = t(int_freq.norm),
                         "Frequentist.profiled" =  t(int_freq.prof),
                         "Bayesian.Quantile" = t(int_bayes_quant),
                         "Bayesian.HPD" = t(int_bayes_hpd)  )
  if(which == 5){
    int_boot_res <- data.frame(ll = boot.ci.Sig_Xi[1, which-4],
                               l = boot.ci.Sig_Xi[1, which-3],
                               m = boot.ci.Sig_Xi[2, which-3],
                               h = boot.ci.Sig_Xi[3, which-4 ],
                               hh = boot.ci.Sig_Xi[3, which-3])
   int_boot_par <- data.frame(ll = boot.ci.Sig_Xi.par[1, which-4],
                             l = boot.ci.Sig_Xi.par[1, which-3],
                             m = boot.ci.Sig_Xi.par[2, which-3],
                             h = boot.ci.Sig_Xi.par[3, which-4 ],
                             hh = boot.ci.Sig_Xi.par[3, which-3])
  }
  else if (which == 4){
    int_boot_res <- data.frame(ll = boot.ci.Sig_Xi[1, which],
                               l = boot.ci.Sig_Xi[1, which-1],
                               m = boot.ci.Sig_Xi[2, which],
                               h = boot.ci.Sig_Xi[3, which ],
                               hh = boot.ci.Sig_Xi[3, which-1])
    int_boot_par <- data.frame(ll = boot.ci.Sig_Xi.par[1, which],
                               l = boot.ci.Sig_Xi.par[1, which-1],
                               m = boot.ci.Sig_Xi.par[2, which],
                               h = boot.ci.Sig_Xi.par[3, which ],
                               hh = boot.ci.Sig_Xi.par[3, which-1])
  }
  if(which == 4 || which == 5) df <- cbind.data.frame(df,
                         "Bootstrap.residual" = t(int_boot_res),
                         "Bootstrap.parametric" = t(int_boot_par))


  mcmc_intervals(df)#, show_density = T) #, rhat = c(1, Rhat[[1]][,2][which]) )
  # Rhat colour by the Rhat (Gelman-Rubin diag) computed above
}

# mu_0
int1 <- intervals_compareByParam(gev_freq_mu0_prof95,
                         gev_freq_mu0_prof75, which = 1) +
  labs(x = "mu0") + theme_piss()
# mu_1
int2 <-intervals_compareByParam(gev_freq_mu1_prof95,
                         gev_freq_mu1_prof75, which = 2) +
  labs(x = "mu1") + theme_piss()
# sigma
int3 <-intervals_compareByParam(gev_freq_sig_prof95,
                         gev_freq_sig_prof75, which = 3) +
  labs(x = "sigma") + theme_piss()
# Xi
int4 <- intervals_compareByParam() +
  labs(x = "xi") + theme_piss(theme_classic()) +
  scale_x_continuous(breaks = c(-0.3, -0.2, -0.1, 0),
                     labels = c(-0.3, -0.2, -0.1, 0))

## All-in-one Sheet
title <- "Confidence intervals comparisons : Frequentists and Bayesians"
grid.arrange(int1, int2, int3, int4, ncol = 2,
             top = grid::textGrob(title,
                                  gp = grid::gpar(col = "#33666C",
                                                  fontsize = 25, hjust = 0.5,
                                                  font = 2), hjust = .5, vjust = 0.4))


####### Posterior Predictive Distribution  ###########

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
n_future <- 116
repl2 <- PissoortThesis::pred_post_samples(n_future = n_future, seed = 12)

post.pred2 <- apply(repl2, 2, function(x) quantile(x, probs = c(0.025,0.5,0.975)))
hpd_pred <- as.data.frame(t(hdi(repl2)))

df.postpred2 <- data.frame(org.data = c(max_years$data,
                                        repl2[sample(10, 1:nrow(repl2)),
                                              117:(116+n_future)] ),
                           q025 = post.pred2["2.5%",], q50 = post.pred2["50%",],
                           q975 = post.pred2["97.5%",], year = 1901:(2016+n_future),
                           'data' = c(rep('original', 116), rep('new', n_future)),
                           hpd.low = hpd_pred$lower, hpd.up = hpd_pred$upper)

col.interval <- c("2.5%-97.5%" = "red", "Median" = "blue2", "HPD 95%" = "green2",
                  "orange", "magenta")
col.data <- c("original" = "cyan", "simulated" = "red", "orange", "magenta")

g.ppd <- ggplot(df.postpred2) +
  geom_line(aes(x = year, y = q025, col = "2.5%-97.5%"), linetype = "dashed") +
  geom_line(aes(x = year, y = q50, col = "Median")) +
  geom_line(aes(x = year, y = q975, col =  "2.5%-97.5%"), linetype = "dashed") +
  geom_line(aes(x = year, y = hpd.low, col = "HPD 95%"), linetype = "dashed") +
  geom_line(aes(x = year, y = hpd.up , col =  "HPD 95%"), linetype = "dashed") +
  geom_vline(xintercept = 2016, linetype = "dashed", size = 0.4, col  = 1) +
  scale_x_continuous(breaks = c(1900, 1950, 2000, 2016, 2050, 2100, 2131),
                     labels = c(1900, 1950, 2000, 2016, 2050, 2100, 2131) ) +
  scale_colour_manual(name = " PP intervals", values = col.interval) +
  geom_point(data = df.postpred2[1:116,],
             aes(x = year, y = org.data), col = "black" ) +
  geom_point(data = df.postpred2[117:nrow(df.postpred2),],
             aes(x = year, y = org.data), col = "orange" ) +
  scale_fill_discrete(name = "Data" ) + #, values = col.data) +
  labs(y = expression( Max~(T~degree*C)), x = "Year",
       title = "Posterior Predictive quantiles with observation + 116 years simulations") +
  theme(legend.position =  c(0.91, 0.12),
        plot.title = element_text(size = 28, colour = "#33666C",
                                  face="bold", hjust = 0.5),
        axis.title = element_text(size = 19, colour = "#33666C", face="bold"),
        legend.title = element_text(size = 19, colour = "#33666C",face="bold") )

## LEngth of the intervals
length.quantil <- df.postpred2$q975 - df.postpred2$q025
length.hpd <- df.postpred2$hpd.up - df.postpred2$hpd.low
df.length.ci <- data.frame(quantiles = length.quantil,
                           hpd = length.hpd,
                           Year = df.postpred2$year)

g.length <- ggplot(df.length.ci) +
  geom_line(aes(x = Year , y = quantiles), col = "red") +
  geom_line(aes(x = Year , y = hpd), col = "green2") +
  labs(title = "Intervals' lengths", y = "Length") +
  scale_x_continuous(breaks = c(1900, 1950, 2000, 2050, 2100, 2131),
                     labels = c(1900, 1950, 2000, 2050, 2100, 2131) ) +
  geom_vline(xintercept = 2016, linetype = "dashed", size = 0.4, col  = 1) +
  theme(plot.title = element_text(size = 17, colour = "#33666C",
                                  face="bold", hjust = 0.5),
        axis.title = element_text(size = 10, colour = "#33666C", face="bold"))

vp <- grid::viewport(width = 0.23,
                    height = 0.28,
                    x = 0.65,
                    y = 0.23)
g.ppd
print(g.length, vp = vp)



# Densities associated with the PPD, with mean(red) and

gg1 <- PissoortThesis::Pred_Dens_ggPlot(1901, repl2)
gg2 <- PissoortThesis::Pred_Dens_ggPlot(1950, repl2)
gg3 <- PissoortThesis::Pred_Dens_ggPlot(2016, repl2)
gg4 <- PissoortThesis::Pred_Dens_ggPlot(2026, repl2)
grid.arrange(gg1, gg2, gg3, gg4, nrow = 1)


## Provide better visualizatons with geom_joy(). ( include it in the package !!!)
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
                           function(x) quantile(x , probs = c(.025, 0.05,
                                                              0.5, 0.95, 0.975)) )
  quantiles_repl2 <- as.data.frame(t(quantiles_repl2))
  quantiles_repl2$year <- colnames(repl2_df)

  repl2_df_gg <- repl2_df[, seq(1, (until + n_future) - from, by = by)] %>%
    gather(year, value)

  col.quantiles <- c("2.5%-97.5%" = "red", "Median" = "black", "HPD 95%" = "green2")

  last_year <- as.character((until + n_future) - from + 1900)
  titl <- paste("Posterior Predictive densities evolution in [ 1901 -", last_year,"] with linear model on location")
  subtitl <- paste("with some quantiles and intervals. The last density is in year ", last_year, ". After 2016 is extrapolation.")

  #browser()
  ## Compute the HPD intervals
  hpd_pred <- as.data.frame(t(hdi(repl2)))
  hpd_pred$year <- colnames(repl2_df)

  g <- ggplot(repl2_df_gg, aes(x = value, y = as.numeric(year) )) +  # %>%rev() inside aes()
    geom_joy(aes(fill = year)) +
    geom_point(aes(x = `2.5%`, y = as.numeric(year), col = "2.5%-97.5%"),
               data = quantiles_repl2, size = 0.9) +
    geom_point(aes(x = `50%`, y = as.numeric(year), col = "Median"),
               data = quantiles_repl2, size = 0.9) +
    geom_point(aes(x = `97.5%`, y = as.numeric(year), col = "2.5%-97.5%"),
               data = quantiles_repl2, size = 0.9) +
    geom_point(aes(x = lower, y = as.numeric(year) , col = "HPD 95%"),
               data = hpd_pred, size = 0.9) +
    geom_point(aes(x = upper, y = as.numeric(year) , col = "HPD 95%"),
               data = hpd_pred, size = 0.9) +
    geom_hline(yintercept = 2016, linetype = "dashed", size = 0.3, col  = 1) +
    scale_fill_viridis(discrete = T, option = "D", direction = -1, begin = .1, end = .9) +
    scale_y_continuous(breaks = c(  seq(1901, 2016, by = by),
      seq(2016, colnames(repl2_df)[ncol(repl2_df)], by = by) ) )  +
    coord_cartesian(xlim = x_coord) +
    theme_piss(theme = theme_minimal()) +
    labs(x = expression( Max~(T~degree*C)), y = "Year",
         title = titl, subtitle = subtitl) +
    scale_colour_manual(name = "Intervals", values = col.quantiles) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(legend.position = c(.952, .37),
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



#####  Return Levels

# Stationary
gibb1$out.chain[, "sigma"] <- exp(gibb1$out.chain[, "logsig"])
rl.post_gg(gibb1$out.chain[, c("mu", "sigma", "xi")], method = "gev")
rl.pred_gg(gibb1$out.chain[, c("mu", "sigma", "xi")], method = "gev",
           qlim = c(30, 40), period = c(1, 5, 15))

#  Nonstationary (linear trend)

# library(evdbayes)
# gibbs.trend$out.chain[,"sigma"] <- exp(gibbs.trend$out.chain[,"logsig"])
# rl.data <-gibbs.trend$out.chain[, c("mu0", "mu1", "sigma", "xi")]
# rl.pred(rl.data, qlim = c(30, 40))
#
# "rl_bayes_trend" <- function(data_year, params, t = 10, m = 10 ){
#   y_m <- -(1 / log(1 - 1/m))
#   t <- seq(max(data_year), max(data_year) + t, 1)
#   rl_m <- (params[1] + params[2] * (t-max(max_years$df$Year))) +
#     (params[3] / params[4]) *  (y_m^params[4] - 1)
#   g <- ggplot(data.frame(r.lvels = rl_m, years = t)) +
#     geom_point(aes(x = years, y = r.lvels))
#   g
#   return(rl_m)
# }
# par_gibbs_trend[3] <- exp(par_gibbs_trend["logsig"])
# rl_bayes_trend(max_years$data,
#                unname(par_gibbs_trend[c("mu0", "mu1", "sigma", "xi")]))

m <- 10   ;   t = 218  ;  start <- 1980
y_m <- -(1 / log(1 - 1/m))
t <- seq(start, start + t, 1)
rl_m <- (tab_quantiles$`50%`[1] + tab_quantiles$`50%`[2] *
           (t-min(max_years$df$Year))) +
  (tab_quantiles$`50%`[5] / tab_quantiles$`50%`[4]) *
  (y_m^tab_quantiles$`50%`[4]- 1)
g <- ggplot(data.frame(r.lvels = rl_m, years = t)) +
  geom_point(aes(x = years, y = r.lvels))
g   ;   rl_m

## Comparing with the frequentist intervals : Retrieve the plot from 2nonstationary.R
gg_rlAll_data + coord_cartesian(ylim = c(30, 40)) +
  geom_point(aes(x = years, y = r.levels),
             data = data.frame(years = t, r.levels = rl_m), col= 2, size = 0.3)




#### Predictive accuracy criterion

## For the Gumbel and the stationary model  computed in
#previous section with Gibbs sampler (gib.mcmc_gumbel  and  gibbs_statio)

# Function for Gumbel model to create mcmc.lists
'mc.listDiag2' <- function (list, subset = c("mu", "logsig")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset]),
            mcmc(list[[4]][, subset])
            )
}

ic_vals.gumb <- gib.mcmc_gumbel$dic.vals
# DIC Values : compute for each chains and then average it.
dic.gum.1 <- mc.listDiag2(gib.mcmc_gumbel$out.ind)[[1]] %>%
    dic_2p(vals = ic_vals.gumb[[1]])
dic.gum.2 <- mc.listDiag2(gib.mcmc_gumbel$out.ind)[[2]] %>%
  dic_2p(vals = ic_vals.gumb[[2]])
dic.gum.3 <- mc.listDiag2(gib.mcmc_gumbel$out.ind)[[3]] %>%
  dic_2p(vals = ic_vals.gumb[[3]])
dic.gum.4 <- mc.listDiag2(gib.mcmc_gumbel$out.ind)[[4]] %>%
  dic_2p(vals = ic_vals.gumb[[4]])
cat(dic.gum.1, dic.gum.2, dic.gum.3, dic.gum.4)

dic.mean.gum <- mean(dic.gum.1, dic.gum.2, dic.gum.3, dic.gum.4)
# WAIC Values
waic.gum.1 <- PissoortThesis::waic( ic_vals.gumb[[1]] )
waic.gum.2 <- PissoortThesis::waic( ic_vals.gumb[[2]] )
waic.gum.3 <- PissoortThesis::waic( ic_vals.gumb[[3]] )
waic.gum.4 <- PissoortThesis::waic( ic_vals.gumb[[4]] )
cat(waic.gum.1, waic.gum.2, waic.gum.3, waic.gum.4)

waic.mean.gum <- mean(waic.gum.1, waic.gum.2, waic.gum.3, waic.gum.4)



# Function for stationary model to create mcmc.lists
'mc.listDiag3' <- function (list, subset = c("mu", "logsig", "xi")) {
  mcmc.list(mcmc(list[[1]][, subset]),
            mcmc(list[[2]][, subset]),
            mcmc(list[[3]][, subset]),
            mcmc(list[[4]][, subset])

  )
}

ic_vals.00 <- gibbs_statio$dic.vals
# DIC Values : compute for each chains and then average it.
dic.statio.1 <- mc.listDiag3(gibbs_statio$out.ind)[[1]] %>%
  PissoortThesis::dic_3p(vals = ic_vals.00[[1]])
dic.statio.2 <- mc.listDiag3(gibbs_statio$out.ind)[[2]] %>%
  PissoortThesis::dic_3p(vals = ic_vals.00[[2]])
dic.statio.3 <- mc.listDiag3(gibbs_statio$out.ind)[[3]] %>%
  PissoortThesis::dic_3p(vals = ic_vals.00[[3]])
dic.statio.4 <- mc.listDiag3(gibbs_statio$out.ind)[[4]] %>%
  PissoortThesis::dic_3p(vals = ic_vals.00[[4]])
cat(dic.statio.1, dic.statio.2, dic.statio.3, dic.statio.4)

dic.mean.statio <- mean(dic.statio.1, dic.statio.2, dic.statio.3, dic.statio.4)
# WAIC Values
waic.statio.1 <- PissoortThesis::waic( ic_vals.00[[1]] )
waic.statio.2 <- PissoortThesis::waic( ic_vals.00[[2]] )
waic.statio.3 <- PissoortThesis::waic( ic_vals.00[[3]] )
waic.statio.4 <- PissoortThesis::waic( ic_vals.00[[4]] )
cat(waic.statio.1, waic.statio.2, waic.statio.3, waic.statio.4)

waic.mean.statio <- mean(waic.statio.1, waic.statio.2, waic.statio.3, waic.statio.4)



## For the nonstationary model computed in this section
ic_vals.trend <- gibbs.trend$dic.vals

# DIC Values : compute for each chains and then average it.
dic.trend.1 <- mc.listDiag4(gibbs.trend$out.ind)[[1]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[1]])
dic.trend.2 <- mc.listDiag4(gibbs.trend$out.ind)[[2]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[2]])
dic.trend.3 <- mc.listDiag4(gibbs.trend$out.ind)[[3]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[3]])
dic.trend.4 <- mc.listDiag4(gibbs.trend$out.ind)[[4]] %>%
  PissoortThesis::dic_4p(vals = ic_vals.trend[[4]])
cat(dic.trend.1, dic.trend.2, dic.trend.3, dic.trend.4)

dic.mean.trend <- mean(dic.trend.1, dic.trend.2, dic.trend.3, dic.trend.4)
# WAIC Values
waic.trend.1 <- PissoortThesis::waic( ic_vals.trend[[1]] )
waic.trend.2 <- PissoortThesis::waic( ic_vals.trend[[2]] )
waic.trend.3 <- PissoortThesis::waic( ic_vals.trend[[3]] )
waic.trend.4 <- PissoortThesis::waic( ic_vals.trend[[4]] )
cat(waic.trend.1, waic.trend.2, waic.trend.3, waic.trend.4)

waic.mean.trend <- mean(waic.trend.1, waic.trend.2, waic.trend.3, waic.trend.4)

## Cross-validation
waic(ic_vals.trend[[1]])




###############################################################################
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
start2 <- list() ;   k <- 1
while(k < 5) { # starting value is randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt2$par, solve(opt2$hessian)))
  svlp <- log_post2(sv[1], sv[2], sv[3], sv[4], sv[5], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start2[[k]] <- sv
    k <- k + 1
  }
}
# k chains with k different starting values
set.seed(101)
gibbs.trend2 <- gibbs.trend2.own(start2, propsd = c(.5, 1.4, 3.5, .2, .15),
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


## Predictive accuracy criterion

ic_vals.trend2 <- gibbs.trend2$dic.vals

# DIC Values. "1:5" takes the 5 parameters of the model, see function
dic.trend2.1 <- mc.listDiag5(gibbs.trend2$out.ind)[[1]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[1]])
dic.trend2.2 <- mc.listDiag5(gibbs.trend2$out.ind)[[2]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[2]])
dic.trend2.3 <- mc.listDiag5(gibbs.trend2$out.ind)[[3]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[3]])
dic.trend2.4 <- mc.listDiag5(gibbs.trend2$out.ind)[[4]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend2[[4]])
cat(dic.trend2.1, dic.trend2.2, dic.trend2.3, dic.trend2.4)

dic.mean.trend2 <- mean(dic.trend2.1, dic.trend2.2, dic.trend2.3,
                        dic.trend2.4)

# WAIC Values
waic.trend2.1 <- PissoortThesis::waic( ic_vals.trend2[[1]] )
waic.trend2.2 <- PissoortThesis::waic( ic_vals.trend2[[2]] )
waic.trend2.3 <- PissoortThesis::waic( ic_vals.trend2[[3]] )
waic.trend2.4 <- PissoortThesis::waic( ic_vals.trend2[[4]] )
cat(waic.trend2.1, waic.trend2.2, waic.trend2.3, waic.trend2.4)
## Interestingly here, all the criterion are lower than for simple trend
# and this models should be preferred... different result from frequentist !!

waic.mean.trend2 <- mean(waic.trend2.1, waic.trend2.2, waic.trend2.3,
                          waic.trend2.4)




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
start3 <- list() ; k <- 1
while(k < 5) { # starting value is randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt3$par, solve(opt3$hessian)))
  svlp <- log_post1(sv[1], sv[2], sv[3], sv[4], sv[5], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start3[[k]] <- sv
    k <- k + 1
  }
}
# k chains with k different starting values
set.seed(101)
gibbs.trend.sig3 <- PissoortThesis::gibbs.trend.sig3own(start3,
                                        propsd = c(.45, 1.3, .2, .5, .1),
                                        iter = 1e3)
colMeans(do.call(rbind, gibbs.trend.sig3$mean_acc.rates))

param.chain3 <- gibbs.trend.sig3$out.chain[, 1:5]


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

dic.sig3.1 <- mc.listDiag4(gibbs.trend.sig3$out.ind, 1:5)[[1]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[1]], sig = T)
dic.sig3.2 <- mc.listDiag4(gibbs.trend.sig3$out.ind, 1:5)[[2]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[2]], sig = T)
dic.sig3.3 <- mc.listDiag4(gibbs.trend.sig3$out.ind, 1:5)[[3]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[3]], sig = T)
dic.sig3.4 <- mc.listDiag4(gibbs.trend.sig3$out.ind, 1:5)[[4]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend.sc[[4]], sig = T)
cat(dic.sig3.1, dic.sig3.2, dic.sig3.3, dic.sig3.4)

dic.mean.sig3 <- mean(dic.sig3.1, dic.sig3.2,
                       dic.sig3.3, dic.sig3.4)
# WAIC Values
waic.sig3.1 <- PissoortThesis::waic( ic_vals.trend.sc[[1]] )
waic.sig3.2 <- PissoortThesis::waic( ic_vals.trend.sc[[2]] )
waic.sig3.3 <- PissoortThesis::waic( ic_vals.trend.sc[[3]] )
waic.sig3.4 <- PissoortThesis::waic( ic_vals.trend.sc[[4]] )
cat(waic.sig3.1, waic.sig3.2, waic.sig3.3, waic.sig3.4)

waic.mean.sig3 <- mean(waic.sig3.1, waic.sig3.2,
                       waic.sig3.3, waic.sig3.4)
# Comparing These values with the ones obtained with simple linear trend, all are
# again for the complex model .....



###############################################################################
## Make the time component cubic in location ?
## mu(t) = mu_0 + mu_1 * t + mu_2 * t^2 + mu_3 * tt^3
##############################################################################


fn33 <- function(par, data) -log_post_mu3(par[1], par[2], par[3],
                                      par[4], par[5], par[6], data )
param33 <- c( mu0 = mean(max_years$df$Max), mu1 = 0, mu2 = 0, mu3 = 0,
             logsig = log(sd(max_years$df$Max)), xi =  -0.1 )
opt33 <- optim(param33, fn33, data = max_years$data,
              method = "BFGS", hessian = T)
opt33

# Starting Values
set.seed(100)
start33 <- list() ;   k <- 1
while(k < 5) { # starting value is randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt33$par, solve(opt33$hessian)))
  svlp <- log_post_mu3(sv[1], sv[2], sv[3], sv[4], sv[5], sv[6], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start33[[k]] <- sv
    k <- k + 1
  }
}

# k chains with k different starting values
set.seed(101)
gibbs.trend33 <- PissoortThesis::gibbs.trend3.own(start33,
                                 propsd = c(.5, 1.4, 3.5, 7, .2, .15),
                                 iter = 1000)
colMeans(do.call(rbind, gibbs.trend33$mean_acc.rates))

param.chain33 <- gibbs.trend33$out.chain[ ,1:5]

## TracePlots
chain.mix <- cbind.data.frame(gibbs.trend33$out.chain,
                              iter.chain = rep(1:500, 4))
mixchains.Own(chain.mix)

## Predictive accuracy criterion

ic_vals.trend33 <- gibbs.trend33$dic.vals

# DIC Values. "1:5" takes the 5 parameters of the model, see function
dic.trend33.1 <- mc.listDiag6(gibbs.trend33$out.ind)[[1]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend33[[1]])
dic.trend33.2 <- mc.listDiag5(gibbs.trend33$out.ind)[[2]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend33[[2]])
dic.trend33.3 <- mc.listDiag5(gibbs.trend33$out.ind)[[3]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend33[[3]])
dic.trend33.4 <- mc.listDiag5(gibbs.trend33$out.ind)[[4]] %>%
  PissoortThesis::dic_5p(vals = ic_vals.trend33[[4]])
cat(dic.trend33.1, dic.trend33.2, dic.trend33.3, dic.trend33.4)

dic.mean.trend33 <- mean(dic.trend33.1, dic.trend33.2, dic.trend33.3,
                        dic.trend33.4)

# WAIC Values
waic.trend33.1 <- PissoortThesis::waic( ic_vals.trend33[[1]] )
waic.trend33.2 <- PissoortThesis::waic( ic_vals.trend33[[2]] )
waic.trend33.3 <- PissoortThesis::waic( ic_vals.trend33[[3]] )
waic.trend33.4 <- PissoortThesis::waic( ic_vals.trend33[[4]] )
cat(waic.trend33.1, waic.trend33.2, waic.trend33.3, waic.trend33.4)
## Interestingly here, all the criterion are lower than for simple trend
# and this models should be preferred... different result from frequentist !!

waic.mean.trend33 <- mean(waic.trend33.1, waic.trend33.2, waic.trend33.3,
                         waic.trend33.4)


##################### Comparisons ("Model Selection")  #######################
#############################################################################

library(loo)

cat(" The following represent the values of the DIC for the models, by ascending complexity",
    dic.mean.gum, dic.mean.statio, dic.mean.trend, dic.mean.trend2,
    dic.mean.sig3, dic.mean.trend33)
cat(" The following represent the values of the WAIC for the same models,",
    waic.mean.gum, waic.mean.statio, waic.mean.trend, waic.mean.trend2,
    waic.mean.sig3, waic.mean.trend33)


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


#save.image("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1_bayes.Rdata")
