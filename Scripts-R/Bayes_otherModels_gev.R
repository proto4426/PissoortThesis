setwd('/home/piss/Documents/Extreme/R resources/IRM')
load("data1.Rdata")

library(PissoortThesis)


######## Make the time component in squared ? ie mu = mu_0 + mu_1 * t^2
##########################################################################

fn2 <- function(par, data) -log_post2(par[1], par[2], par[3], par[4], par[5], data )
param2 <- c( mu0 = mean(max_years$df$Max), mu1 = 0, mu2 = 0,
            logsig = log(sd(max_years$df$Max)), xi =  -0.1 )
opt2 <- optim(param2, fn2, data = max_years$data,
             method = "BFGS", hessian = T)
opt2

# Starting Values
set.seed(100)
start <- list() ; k <- 1
while(k < 5) { # starting value is randomly selected from a distribution
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
set.seed(100)
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



## Predictive accuracy criterion

ic_vals <- gibbs.trend2$dic.vals
'dic2' <- function(out, vals) {
  pm <- colMeans(out) ;   pmv <- log_post2(pm[1], pm[2], pm[3], pm[4], pm[5], data)
  pmv <- sum(pmv, na.rm = TRUE) ;   vec1 <- rowSums(vals, na.rm = TRUE)
  2*pmv - 4*mean(vec1)
}
# DIC Values. "1:5" takes the 5 parameters of the model, see function
dic2( mc.listDiag(gibbs.trend2$out.ind, 1:5)[[1]],
     ic_vals[[1]] ) - 2200
dic2( mc.listDiag(gibbs.trend2$out.ind, 1:5)[[2]],
     ic_vals[[1]] ) - 2200
dic2( mc.listDiag(gibbs.trend2$out.ind, 1:5)[[3]],
     ic_vals[[1]] ) - 2200
dic2( mc.listDiag(gibbs.trend2$out.ind, 1:5)[[4]],
     ic_vals[[1]] ) - 2200
# WAIC Values
waic( ic_vals[[1]] ) - 2200
waic( ic_vals[[2]]) - 2200
waic( ic_vals[[3]]) - 2200
waic( ic_vals[[4]] ) - 2200
## Interestingly here, all the criterion are lower than for simple trend
# and this models should be preferred... different result from frequentist !!


## Cross-validation





######## Model with inear trend + varying scale parameter (exp link)
########   mu = mu0 + mu1 * t,     logsig = sig0 + sig1 * tt
##########################################################################

fn3 <- function(par, data) -log_post3(par[1], par[2], par[3], par[4], par[5], data )
param3 <- c( mu0 = mean(max_years$df$Max), mu1 = 0,
             sig0 = log(sd(max_years$df$Max)), sig1 = 0,  xi =  -0.1 )
opt3 <- optim(param3, fn3, data = max_years$data,
              method = "BFGS", hessian = T)
opt3

# Starting Values
set.seed(100)
start <- list() ; k <- 1
while(k < 5) { # starting value is randomly selected from a distribution
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
set.seed(100)
gibbs.trend.sig3 <- gibbs.trend.sig3own(start, propsd = c(.45, 1.3, .2, .5, .1),
                                     iter = 1000)
colMeans(do.call(rbind, gibbs.trend.sig3$mean_acc.rates))

param.chain3 <- gibbs.trend.sig3$out.chain[ ,1:5]


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

ic_vals <- gibbs.trend.sig3$dic.vals
'dic3' <- function(out, vals) {
  pm <- colMeans(out) ;   pmv <- log_post3(pm[1], pm[2], pm[3], pm[4], pm[5], data)
  pmv <- sum(pmv, na.rm = TRUE) ;   vec1 <- rowSums(vals, na.rm = TRUE)
  2*pmv - 4*mean(vec1)
}
# DIC Values. "1:5" takes the 5 parameters of the model, see function
dic3( mc.listDiag(gibbs.trend.sig3$out.ind, 1:5)[[1]], ic_vals[[1]] )
dic3( mc.listDiag(gibbs.trend.sig3$out.ind, 1:5)[[2]], ic_vals[[2]] )
dic3( mc.listDiag(gibbs.trend.sig3$out.ind, 1:5)[[3]], ic_vals[[3]] )
dic3( mc.listDiag(gibbs.trend.sig3$out.ind, 1:5)[[4]], ic_vals[[4]] )
# WAIC Values
waic( ic_vals[[1]] )
waic( ic_vals[[2]])
waic( ic_vals[[3]])
waic( ic_vals[[4]] )

# Comparing These values with the ones obtained with simple linear trend, all are
# again for the complex model .....


## Cross-validation

library(loo)

extract_log_lik()

