
# ===============================================================
#' @title Negative GEV log-likelihood and Log-Posterior with
#' uninformative normal priors.
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Functions to compute the negative log-likelihood of a GEV
#' distribution. Whereas this can be used for other purpose, this
#' also allows in particular to compute the log-posterior with diffuse normal
#' priors. Note that we could add parameters to control the informativness of
#' the priors, but as we have no reliable information, we decide to arbitrarily
#'  fix it to large values, to improve computation . (More parameters are
#'  harmful for computation time)
#'
#' @param mu numeric representig the location parameter of the GEV
#' @param sig or \code{logsig} are numeric representig the scale parameter
#' of the GEV. BE CAREFULL : this is in logarithm for the log_post0 only in
#' order to avoid computational problems later in the MCMC's.
#' @param xi numeric representig the shape parameter of the GEV
#' @param data numeric vector representing the data (GEV) of interest.
#'
#' @return a numeric value representing the negative log-likelihood or the log-posterior
#' of interest.
#' @examples
#' data('max_years')
#' # Optimize the log-Posterior Density Function to find starting values
#' fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
#' param <- c(mean(max_years$df$Max),log(sd(max_years$df$Max)), 0.1 )
#' opt <- nlm(fn, param, data = max_years$data,
#'            hessian=T, iterlim = 1e5)
#' start <- opt$estimate
#'
#' @rdname log_post0
#' @export
'gev.nloglik' = function(mu, sig, xi, data){
  y = 1 + (xi * (data - mu))/sig
  if((sig < 0) || (min(y) < 0) || (is.na(y))) {
    ans = 1e+06
  } else {
    term1 = length(data) * logb(sig)
    term2 = sum((1 + 1/xi) * logb(y))
    term3 = sum(y^(-1/xi))
    ans = term1 + term2 + term3
  }
  ans
}
#' @rdname log_post0
#' @export
'log_post0' <- function(mu, logsig,xi, data) {
  # Posterior Density Function
  # Compute the log_posterior in a stationary context.
  # Be careful to incorporate the fact that the distribution can have finite endpoints.
  llhd <- -(gev.nloglik(mu = mu, sig = exp(logsig),
                        xi = xi, data = data))
  lprior <- dnorm(mu, sd = 50, log = TRUE)
  lprior <- lprior + dnorm(logsig, sd = 50, log = TRUE)
  lprior <- lprior + dnorm(xi, sd = 5, log = TRUE)
  lprior + llhd
}




# ===============================================================
#' @title Plots to assess the mixing of the Chains
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute the ggplots for each parameter of interest in a single page.
#'
#' @param data numeric vector containing the GEV block-maxima
#' @param ... Other parameters from \code{gridExtra::grid.arrange()}
#' @param title Global title for the plot
#' @return a grid.arrange() of ggplots.
#' @examples
#' data("max_years")
#' fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
#' param <- c(mean(max_years$df$Max),log(sd(max_years$df$Max)), 0.1 )
#' # opt <- optim(param, fn, data = max_years$data,
#' #              method="BFGS", hessian = TRUE)
#' opt <- nlm(fn, param, data = max_years$data,
#'            hessian=T, iterlim = 1e5)
#' start <- opt$estimate
#' Sig <- solve(opt$hessian)
#' ev <- eigen( (2.4/sqrt(2))^2 * Sig)
#' varmat <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
#' # (MH)
#' set.seed(100)
#' mh.mcmc1 <- MH_mcmc.own(start, varmat %*% c(.1,.3,.4))
#' mh.mcmc1$mean.acc_rates
#'
#' chains.plotOwn(mh.mcmc1$out.chain)
#'
#' # (GIBBS)
#' # k chains with k different starting values
#' set.seed(100)
#' gibbs.trend <- gibbs.trend.own(start, propsd = c(.5, 1.9, .15, .12),
#'                               iter = 1000)
#'## TracePlots
#' chain.mix <- cbind.data.frame(gibbs.trend$out.chain,
#'                              iter.chain = rep(1:500, 4))
#' mixchains.Own(chain.mix)
#' @import ggplot2  gridExtra grid
#' @rdname ggplotbayesfuns
#' @export
'chains.plotOwn' <- function(data, ...,
                             title = "TracePlots of the generated Chains " ){
  grid.arrange(
    ggplot(data) + geom_line(aes(x = iter, y = mu)) + theme_piss(18,16) +
      labs(ylab = expression(mu)),
    ggplot(data) + geom_line(aes(x = iter, y = logsig)) + theme_piss(18,16),
    ggplot(data) + geom_line(aes(x = iter, y = xi)) + theme_piss(18,16), ... ,
    ncol = 1,
    top = textGrob(title,
                   gp = gpar(col ="darkolivegreen4",
                             fontsize = 25, font = 4))
  )
}
#' @rdname ggplotbayesfuns
#' @export
'mixchains.Own' <- function(data, moreplot = F,
                             title = "TracePlots of the generated Chains " ){
  grid.arrange(
    ggplot(data, aes(x = iter.chain, y = mu, col = as.factor(chain.nbr))) +
      geom_line() + theme_piss(18,16, theme = theme_classic()) +
      scale_colour_brewer(name = "chain nr", palette = "Set1") +
      guides(colour = guide_legend(override.aes = list(size= 1.2))),
    ggplot(data, aes(x = iter.chain, y = logsig, col = as.factor(chain.nbr))) +
      geom_line() + theme_piss(18,16, theme = theme_classic()) +
      scale_colour_brewer(name = "chain nr", palette = "Set1") +
      guides(colour = guide_legend(override.aes = list(size= 1.2))),
    ggplot(data, aes(x = iter.chain, y = xi, col = as.factor(chain.nbr))) +
      geom_line() + theme_piss(18,16, theme = theme_classic()) +
      scale_colour_brewer(name = "chain nr", palette = "Set1") +
      guides(colour = guide_legend(override.aes = list(size= 1.2))),
    ncol = 1, top = textGrob(title,
                   gp = gpar(col ="darkolivegreen4",
                             fontsize = 25, font = 4))
  )
}

# ===============================================================
#' @export MH_mcmc.own
#' @title Metropolis-Hastings algorithm for GEV
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute return levels plot of nonstationary model with the data (in years)
#'
#' @param start numeric vector of length 3 containing the starting values for the parameters theta=
#'(location, LOG-scale and shape). It is advised explore different ones, and typically take the MPLE
#' @param varmat.prop The proposal's variance : controlling the cceptance rate.
#' To facilitate convergence, it
#' is advised to target an acceptance rate of around 0.25 when all components of theta are updated
#' simultaneously, and 0.40 when the components are updated one at a time.
#' @param data  numeric vector containing the GEV in block-maxima
#' @param iter The number of iterations of the algorithm. Must e high enough to ensure convergence
#'
#' @return A named list containing
#' \describe{
#' \item{\code{mean.acc_rates} : the mean of the acceptance rates}
#' \item{\code{out.chain} : The generated chain}
#' }
#' @examples
#' data("max_years")
#' fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
#' param <- c(mean(max_years$df$Max),log(sd(max_years$df$Max)), 0.1 )
#' # opt <- optim(param, fn, data = max_years$data,
#' #              method="BFGS", hessian = TRUE)
#' opt <- nlm(fn, param, data = max_years$data,
#'            hessian=T, iterlim = 1e5)
#' start <- opt$estimate
#' Sig <- solve(opt$hessian)
#' ev <- eigen( (2.4/sqrt(2))^2 * Sig)
#' varmat <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
#' # (MH)
#' set.seed(100)
#' mh.mcmc1 <- MH_mcmc.own(start, varmat %*% c(.1,.3,.4))
'MH_mcmc.own' <- function(start, varmat.prop,
                          data = max_years$data, iter = 2000){

  out <- matrix(NA, nrow = iter+1, ncol = 3)
  dimnames(out) <- list(1:(iter+1), c("mu", "logsig", "xi"))
  out[1,] <- start
  lpost_old <- log_post0(out[1,1], out[1,2], out[1,3], data)
  if(!is.finite(lpost_old))
    stop("starting values give non-finite log_post")
  acc_rates <- numeric(iter)
  for(t in 1:iter) {
    prop <- rnorm(3) %*% varmat.prop + out[t,]  # add tuning parameter delta ?
    lpost_prop <- log_post0(prop[1], prop[2], prop[3], data)
    r <- exp(lpost_prop - lpost_old) # as the prop is symmetric
    if(r > runif(1)) {
      out[t+1,] <- prop
      lpost_old <- lpost_prop
    }
    else out[t+1,] <- out[t,]
    acc_rates[t] <- min(r, 1)
  }
  return(list(mean.acc_rates = mean(acc_rates),
              out.chain = data.frame(out, iter = 1:(iter+1))))
}


# ===============================================================
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @title Gibbs Sampler for GEV (MCMC)
#' @description
#' Compute return levels plot of nonstationary model with the data (in years)
#'
#' @param start numeric vector of length 3 containing the starting values for the parameters theta=
#'(location, LOG-scale and shape). It is advised explore different ones, and typically take the MPLE
#' @param proposd The proposal's standard deviations : controlling the cceptance rate.
#' To facilitate convergence, it is advised to target an acceptance rate of around 0.25
#' when all components of theta are updated  simultaneously,
#' and 0.40 when the components are updated one at a time.
#' is as proposed by Coles (2001) but we should tune this value.
#'  (from package \code{ismev})
#' @param data  numeric vector containing the GEV in block-maxima
#' @param iter The number of iterations of the algorithm. Must e high enough to ensure convergence
#'
#' @return A named list containing
#' \describe{
#' \item{\code{mean.acc_rates} : the meanS of the acceptance rates}
#' \item{\code{out.chain} : The generated chainS}
#' }
#' @examples
#' data("max_years")
#' fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
#' param <- c(mean(max_years$df$Max),log(sd(max_years$df$Max)), 0.1 )
#' # opt <- optim(param, fn, data = max_years$data,
#' #              method="BFGS", hessian = TRUE)
#' opt <- nlm(fn, param, data = max_years$data,
#'            hessian=T, iterlim = 1e5)
#' start <- opt$estimate
#' Sig <- solve(opt$hessian)
#' ev <- eigen( (2.4/sqrt(2))^2 * Sig)
#' varmat <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
#' ## GIBBS
#'  set.seed(100)
#' iter <- 2000
#' gibb1 <- gibbs_mcmc.own(start, iter = iter)
#' @export
"gibbs_mcmc.own" <- function (start , propsd = c(.4, .1, .1),
                              iter = 2000, data = max_years$data ) {

 out <- data.frame(mu = rep(NA, iter+1),
                  logsig = rep(NA, iter+1),
                  xi = rep(NA, iter+1))

 out[1,] <- start
 out <- cbind.data.frame(out, iter = 1:(iter+1))
 lpost_old <- log_post0(out[1,1], out[1,2], out[1,3], data)
 if(!is.finite(lpost_old))
   stop("starting values give non-finite log_post")
 acc_rates <- matrix(NA, nrow = iter, ncol = 3)

 data <- max_years$data
 for (t in 1:iter) {
   prop1 <- rnorm(1, mean = out[t,1], propsd[1]) # symmetric too
   # so that it removes in the ratio.

   lpost_prop <- log_post0(prop1, out[t,2], out[t,3], data)
   r <- exp(lpost_prop - lpost_old)
   if(r > runif(1)) {
     out[t+1,1] <- prop1
     lpost_old <- lpost_prop
   }
   else out[t+1,1] <- out[t,1]
   acc_rates[t,1] <- min(r, 1)

   prop2 <- rnorm(1, mean = out[t,2], propsd[2])
   lpost_prop <- log_post0(out[t+1,1], prop2, out[t,3], data)
   r <- exp(lpost_prop - lpost_old)
   if(r > runif(1)) {
     out[t+1,2] <- prop2
     lpost_old <- lpost_prop
   }
   else out[t+1,2] <- out[t,2]
   acc_rates[t,2] <- min(r, 1)

   prop3 <- rnorm(1, mean = out[t,3], propsd[3])
   lpost_prop <- log_post0(out[t+1,1],out[t+1,2], prop3, data)
   r <- exp(lpost_prop - lpost_old)
   if(r > runif(1)) {
     out[t+1,3] <- prop3
     lpost_old <- lpost_prop
   }
   else out[t+1,3] <- out[t,3]
   acc_rates[t,3] <- min(r, 1)
 }
 return(list(mean.acc_rates = apply(acc_rates, 2, mean),
             out.chain = out))
}




# ===============================================================
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @title Gibbs Sampler for nonstationary GEV (MCMC)
#' @description
#' Compute the Gibbs sampler accounting for nonstationarity (trend) in GEV. It is computed
#' based on a diffuse normal prior.
#' @param start named list of length 4 (this number determines the number of chains generated)
#' containing the starting values for the parameters theta=(intercept mu0, trend mu1,
#' LOG-scale and shape).
#' It is advised explore different ones, and typically take the MPLE.
#' @param proposd The proposal's standard deviations : controlling the cceptance rate.
#' To facilitate convergence, it is advised to target an acceptance rate of around 0.25
#' when all components of theta are updated  simultaneously,
#' and 0.40 when the components are updated one at a time.
#' It must be wel chosen (e.g. Trial-and-error method)
#' @param data  numeric vector containing the GEV in block-maxima
#' @param iter The number of iterations of the algorithm. Must e high enough to ensure convergence
#'
#' @return A named list containing
#' \describe{
#' \item{\code{n.chains} : The number of chains generated melted in a data.frame}
#' \item{\code{mean.acc_rates} : the meanS of the acceptance rates}
#' \item{\code{out.chain} : The generated chainS}
#' \item{\code{dic.vals} : contains the DIC values (for further diagnostics on
#' predictive accuracy, see ?dic)}
#' \item{\code{out.ind} : The generated individual chainS (in a list)}

#' }
#' @examples
#' data("max_years")
#' data <- max_years$data
#'
#' fn <- function(par, data) -log_post1(par[1], par[2], par[3],
#'                                      par[4], data)
#' param <- c(mean(max_years$df$Max), 0, log(sd(max_years$df$Max)), -0.1 )
#' opt <- optim(param, fn, data = max_years$data,
#'              method = "BFGS", hessian = T)
#'
#' # Starting Values
#' set.seed(100)
#' start <- list() ; k <- 1
#' while(k < 5) { # starting value is randomly selected from a distribution
#'   # that is overdispersed relative to the target
#'   sv <- as.numeric(rmvnorm(1, opt$par, 50 * solve(opt$hessian)))
#'   svlp <- log_post1(sv[1], sv[2], sv[3], sv[4], max_years$data)
#'   print(svlp)
#'   if(is.finite(svlp)) {
#'     start[[k]] <- sv
#'     k <- k + 1
#'   }
#' }
#'
#' # k chains with k different starting values
#' set.seed(100)
#' gibbs.trend <- gibbs.trend.own(start, propsd = c(.5, 1.9, .15, .12),
#'                                iter = 1000)
#' colMeans(do.call(rbind, gibbs.trend$mean_acc.rates))
#'
#' param.chain <- gibbs.trend$out.chain[ ,1:4]
#'
#' ### Plot of the chains
#' chains.plotOwn(gibbs.trend$out.chain )
#' @export
'log_post1' <- function(mu0, mu1, logsig, xi, data,
                        model.mu = mu0 + mu1 * tt,
                        mnpr = c(30,0,0,0), sdpr = c(40,40,10,10)) {
  theta <- c(mu0, mu1, logsig, xi)
  tt <- ( min(max_years$df$Year):max(max_years$df$Year) -
             mean(max_years$df$Year) ) / length(max_years$data)
  mu <-  model.mu
  llhd1 <- sum(evd::dgev(data, loc = mu, scale = exp(logsig), xi,
                         log = TRUE))
  lprior <- sum(dnorm(theta, mean = mnpr, sd = sdpr, log = TRUE))
  lprior + llhd1 #+ llhd2
}
#' @export
'gibbs.trend.own' <- function (start, propsd = c(.5, 2.5, .08, .08),
                               iter = 1000, data = max_years$data) {
  # To store values inside
  acc_rate.list <- list() ;  ic_val.list <- list() ;  out.ind <- list()

  hf <- ceiling(iter/2 + 1) # Determines values for burn.in (see end)

  out.fin <- data.frame(mu = numeric(0),
                        mu1 = numeric(0),
                        logsig = numeric(0),
                        xi = numeric(0),
                        chain.nbr = character(0))
  nr.chain <- length(start)   ;    time <- proc.time()

 for(k in 1:nr.chain) {
   out <- data.frame(mu = rep(NA, iter+1),
                     mu1 = rep(NA, iter+1),
                     logsig = rep(NA, iter+1),
                     xi = rep(NA, iter+1))

   out[1,] <- start[[k]]
   out <- cbind.data.frame(out, chain.nbr = rep(as.factor(k), iter+1))

   lpost_old <- log_post1(out[1,1], out[1,2], out[1,3], out[1,4], data)

   # For DIC computation
   ic_vals <- matrix(NA, nrow = iter+1, ncol = length(data))
   ic_vals[1,] <- log_post1(out[1,1], out[1, 2], out[1,3], out[1,4],
                            data)

   if(!is.finite(lpost_old))
     stop("starting values give non-finite log_post")
    acc_rates <- matrix(NA, nrow = iter, ncol = 4)

   for(t in 1:iter) {
     prop1 <- rnorm(1, mean = out[t,1], propsd[1])
     lpost_prop <- log_post1(prop1, out[t,2], out[t,3], out[t,4], data)
     r <- exp(lpost_prop - lpost_old)
     if(r > runif(1)) {
       out[t+1,1] <- prop1
       lpost_old <- lpost_prop
     }
     else out[t+1,1] <- out[t,1]
     acc_rates[t,1] <- min(r, 1)

     prop2 <- rnorm(1, mean = out[t,2], propsd[2])
     lpost_prop <- log_post1(out[t+1,1], prop2, out[t,3], out[t,4], data)
     r <- exp(lpost_prop - lpost_old)
     if(r > runif(1)) {
       out[t+1,2] <- prop2
       lpost_old <- lpost_prop
     }
     else out[t+1,2] <- out[t,2]
     acc_rates[t,2] <- min(r, 1)

     prop3 <- rnorm(1, mean = out[t,3], propsd[3])
     lpost_prop <- log_post1(out[t+1,1], out[t+1,2], prop3, out[t,4], data)
     r <- exp(lpost_prop - lpost_old)
     if(r > runif(1)) {
       out[t+1,3] <- prop3
       lpost_old <- lpost_prop
     }
     else out[t+1,3] <- out[t,3]
     acc_rates[t,3] <- min(r, 1)

     prop4 <- rnorm(1, mean = out[t,4], propsd[4])
     lpost_prop <- log_post1(out[t+1,1], out[t+1,2], out[t+1,3], prop4, data)
     r <- exp(lpost_prop - lpost_old)
     if(r > runif(1)) {
       out[t+1,4] <- prop4
       lpost_old <- lpost_prop
     }
     else out[t+1,4] <- out[t,4]
     acc_rates[t,4] <- min(r, 1)

     # For DIC
     ic_vals[t+1,] <- log_post1(out[1,1], out[1, 2], out[1,3], out[1,4],
                                data)
   }
   acc_rate.list[[k]] <- apply(acc_rates, 2, mean )
   ic_val.list[[k]] <- ic_vals[-(1:hf), ]
   out.ind[[k]] <- out

   # Combine Chains And Remove Burn-In Period
   out.fin <- rbind.data.frame(out.fin, out[-(1:hf), ])
  # out.fin <- cbind.data.frame(out.fin)
                              # chain.nmbr = rep(k, nrow(out.fin)))

   print(paste("time is ", round((proc.time() - time)[3], 5), " sec"))
 }

  out <- cbind.data.frame(out.fin,
                          iter = (1:nrow(out.fin)))

 return( list(n.chains = length(start),
              mean_acc.rates = acc_rate.list,
              out.chain = out,
              dic.vals = ic_val.list,
              out.ind = out.ind) )
}



# ===============================================================
#' @name gibbs_trend2
#' @title Return Levels with nonstationarity
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute return levels plot of nonstationary model with the data (in years)
#'
#' @param data numeric vector containing the GEV block-maxima
#' @param gev_nstatio Nonstationary GEV fitted model of class \code{gev.fit}
#' (from package \code{ismev})
#' @param t Maximum time period for which the return levels are considered
#' @param m Return period
#'
#' @return The return levels for the considered time period (t)
#' @examples
#' rl_10_lin <- return.lvl.nstatio(max_years$df$Year,
#' @rdname gibbs2
#' @export
'log_post2' <- function(mu0, mu1, mu2, logsig, xi, data,
                        model.mu = mu0 + mu1 * tt + mu2 * tt^2,
                        mnpr = c(30,0,0,0,0), sdpr = c(40,40,10,10,10)) {
  theta <- c(mu0, mu1, mu2, logsig, xi)
  tt <- ( min(max_years$df$Year):max(max_years$df$Year) -
            mean(max_years$df$Year) ) / length(max_years$data)
  mu <-  model.mu
  llhd1 <- sum(evd::dgev(data, loc = mu, scale = exp(logsig), xi,
                         log = TRUE))
  lprior <- sum(dnorm(theta, mean = mnpr, sd = sdpr, log = TRUE))
  lprior + llhd1 #+ llhd2
}
#' @rdname gibbs2
#' @export
'gibbs.trend2.own' <- function (start, propsd = c(.5, 2.5, 2, .08, .08),
                               iter = 1000, data = max_years$data) {
  # To store values inside
  acc_rate.list <- list() ;  ic_val.list <- list() ;  out.ind <- list()

  hf <- ceiling(iter/2 + 1) # Determines values for burn.in (see end)

  out.fin <- data.frame(mu = numeric(0),
                        mu1 = numeric(0),
                        mu2 = numeric(0),
                        logsig = numeric(0),
                        xi = numeric(0),
                        chain.nbr = character(0))
  nr.chain <- length(start)   ;    time <- proc.time() ;  k = 1
  while(k <= nr.chain) {
    out <- data.frame(mu = rep(NA, iter+1),
                      mu1 = rep(NA, iter+1),
                      mu2 = rep(NA, iter+1),
                      logsig = rep(NA, iter+1),
                      xi = rep(NA, iter+1))

    out[1,] <- start[[k]]
    out <- cbind.data.frame(out, chain.nbr = rep(as.factor(k), iter+1))

    lpost_old <- log_post2(out[1,1], out[1,2], out[1,3], out[1,4], out[1,5], data)

    # For DIC computation
    ic_vals <- matrix(NA, nrow = iter+1, ncol = length(data))
    ic_vals[1,] <- log_post2(out[1,1], out[1, 2], out[1,3], out[1,4],  out[1,5],
                             data)

    if(!is.finite(lpost_old))
      stop("starting values give non-finite log_post")
    acc_rates <- matrix(NA, nrow = iter, ncol = 5)

    for(t in 1:iter) {
      prop1 <- rnorm(1, mean = out[t,1], propsd[1])
      lpost_prop <- log_post2(prop1, out[t,2], out[t,3], out[t,4],  out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,1] <- prop1
        lpost_old <- lpost_prop
      }
      else out[t+1,1] <- out[t,1]
      acc_rates[t,1] <- min(r, 1)

      prop2 <- rnorm(1, mean = out[t,2], propsd[2])
      lpost_prop <- log_post2(out[t+1,1], prop2, out[t,3], out[t,4],  out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,2] <- prop2
        lpost_old <- lpost_prop
      }
      else out[t+1,2] <- out[t,2]
      acc_rates[t,2] <- min(r, 1)

      prop3 <- rnorm(1, mean = out[t,3], propsd[3])
      lpost_prop <- log_post2(out[t+1,1], out[t+1,2], prop3, out[t,4],
                              out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,3] <- prop3
        lpost_old <- lpost_prop
      }
      else out[t+1,3] <- out[t,3]
      acc_rates[t,3] <- min(r, 1)

      prop4 <- rnorm(1, mean = out[t,4], propsd[4])
      lpost_prop <- log_post2(out[t+1,1], out[t+1,2], out[t+1,3], prop4,
                              out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,4] <- prop4
        lpost_old <- lpost_prop
      }
      else out[t+1,4] <- out[t,4]
      acc_rates[t,4] <- min(r, 1)

      prop5 <- rnorm(1, mean = out[t,5], propsd[5])
      lpost_prop <- log_post2(out[t+1,1], out[t+1,2], out[t+1,3], out[t+1,4],
                              prop5, data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,5] <- prop5
        lpost_old <- lpost_prop
      }
      else out[t+1,5] <- out[t,5]
      acc_rates[t,5] <- min(r, 1)

      # For DIC
      ic_vals[t+1,] <- log_post2(out[1,1], out[1, 2], out[1,3], out[1,4], out[1,5],
                                 data)
    }
    acc_rate.list[[k]] <- apply(acc_rates, 2, mean )
    ic_val.list[[k]] <- ic_vals[-(1:hf), ]
    out.ind[[k]] <- out

    # Combine Chains And Remove Burn-In Period
    out.fin <- rbind.data.frame(out.fin, out[-(1:hf), ])
    # out.fin <- cbind.data.frame(out.fin)
    # chain.nmbr = rep(k, nrow(out.fin)))

    print(paste("time is ",round((proc.time() - time)[3], 5), " sec"))
    k = k +1
  }

  out <- cbind.data.frame(out.fin,
                          iter = (1:nrow(out.fin)))

  return( list(n.chains = length(start),
               mean_acc.rates = acc_rate.list,
               out.chain = out,
               dic.vals = ic_val.list,
               out.ind = out.ind) )
}


# ===============================================================
#' @name gibbstrend3
#' @title Return Levels with nonstationarity
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute return levels plot of nonstationary model with the data (in years)
#'
#' @param data numeric vector containing the GEV block-maxima
#' @param gev_nstatio Nonstationary GEV fitted model of class \code{gev.fit}
#' (from package \code{ismev})
#' @param t Maximum time period for which the return levels are considered
#' @param m Return period
#'
#' @return The return levels for the considered time period (t)
#' @examples
#' rl_10_lin <- return.lvl.nstatio(max_years$df$Year,
#'

## for the model 3 :allowing linear trend in mu and log-linked varying scale param.
#' @rdname gibbs3
#' @export
'log_post3' <- function(mu0, mu1, sig0, sig1, xi, data,
                        model.mu = mu0 + mu1 * tt,
                        mnpr = c(30,0,1,0, 0), sdpr = c(10,40,10,10, 10)) {
  theta <- c(mu0, mu1, sig0, sig1, xi)
  tt <- ( min(max_years$df$Year):max(max_years$df$Year) -
            mean(max_years$df$Year) ) / length(max_years$data)
  mu <-  model.mu
  logsig <- sig0 + sig1 * tt
  llhd1 <- sum(evd::dgev(data, loc = mu, scale = exp(logsig), xi,
                         log = TRUE))
  lprior <- sum(dnorm(theta, mean = mnpr, sd = sdpr, log = TRUE))
  lprior + llhd1 #+ llhd2
}
#' @rdname gibbs3
#' @export
'gibbs.trend.sig3own' <- function (start, propsd = c(.5, 2.5, 2, .08, .08),
                                iter = 1000, data = max_years$data) {
  # To store values inside
  acc_rate.list <- list() ;  ic_val.list <- list() ;  out.ind <- list()

  hf <- ceiling(iter/2 + 1) # Determines values for burn.in (see end)

  out.fin <- data.frame(mu = numeric(0),
                        mu1 = numeric(0),
                        sig0 = numeric(0),
                        sig1 = numeric(0),
                        xi = numeric(0),
                        chain.nbr = character(0))
  nr.chain <- length(start)   ;    time <- proc.time()  ;  k = 1
  while(k <= nr.chain) {
    out <- data.frame(mu = rep(NA, iter+1),
                      mu1 = rep(NA, iter+1),
                      sig0 = rep(NA, iter+1),
                      sig1 = rep(NA, iter+1),
                      xi = rep(NA, iter+1))

    out[1,] <- start[[k]]
    out <- cbind.data.frame(out, chain.nbr = rep(as.factor(k), iter+1))

    lpost_old <- log_post3(out[1,1], out[1,2], out[1,3], out[1,4], out[1,5], data)

    # For DIC computation
    ic_vals <- matrix(NA, nrow = iter+1, ncol = length(data))
    ic_vals[1,] <- log_post3(out[1,1], out[1, 2], out[1,3], out[1,4],  out[1,5],
                             data)

    if(!is.finite(lpost_old))
      stop("starting values give non-finite log_post")
    acc_rates <- matrix(NA, nrow = iter, ncol = 5)

    for(t in 1:iter) {
      prop1 <- rnorm(1, mean = out[t,1], propsd[1])
      lpost_prop <- log_post3(prop1, out[t,2], out[t,3], out[t,4],  out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,1] <- prop1
        lpost_old <- lpost_prop
      }
      else out[t+1,1] <- out[t,1]
      acc_rates[t,1] <- min(r, 1)

      prop2 <- rnorm(1, mean = out[t,2], propsd[2])
      lpost_prop <- log_post3(out[t+1,1], prop2, out[t,3], out[t,4],  out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,2] <- prop2
        lpost_old <- lpost_prop
      }
      else out[t+1,2] <- out[t,2]
      acc_rates[t,2] <- min(r, 1)

      prop3 <- rnorm(1, mean = out[t,3], propsd[3])
      lpost_prop <- log_post3(out[t+1,1], out[t+1,2], prop3, out[t,4],
                              out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,3] <- prop3
        lpost_old <- lpost_prop
      }
      else out[t+1,3] <- out[t,3]
      acc_rates[t,3] <- min(r, 1)

      prop4 <- rnorm(1, mean = out[t,4], propsd[4])
      lpost_prop <- log_post3(out[t+1,1], out[t+1,2], out[t+1,3], prop4,
                              out[t,5], data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,4] <- prop4
        lpost_old <- lpost_prop
      }
      else out[t+1,4] <- out[t,4]
      acc_rates[t,4] <- min(r, 1)

      prop5 <- rnorm(1, mean = out[t,5], propsd[5])
      lpost_prop <- log_post3(out[t+1,1], out[t+1,2], out[t+1,3], out[t+1,4],
                              prop5, data)
      r <- exp(lpost_prop - lpost_old)
      if(r > runif(1)) {
        out[t+1,5] <- prop5
        lpost_old <- lpost_prop
      }
      else out[t+1,5] <- out[t,5]
      acc_rates[t,5] <- min(r, 1)

      # For DIC
      ic_vals[t+1,] <- log_post3(out[1,1], out[1, 2], out[1,3], out[1,4], out[1,5],
                                 data)
    }
    acc_rate.list[[k]] <- apply(acc_rates, 2, mean )
    ic_val.list[[k]] <- ic_vals[-(1:hf), ]
    out.ind[[k]] <- out

    # Combine Chains And Remove Burn-In Period
    out.fin <- rbind.data.frame(out.fin, out[-(1:hf), ])
    # out.fin <- cbind.data.frame(out.fin)
    # chain.nmbr = rep(k, nrow(out.fin)))

    print(paste("time is ",round((proc.time() - time)[3], 5), " sec"))
    k = k + 1
  }

  out <- cbind.data.frame(out.fin,
                          iter = (1:nrow(out.fin)))

  return( list(n.chains = length(start),
               mean_acc.rates = acc_rate.list,
               out.chain = out,
               dic.vals = ic_val.list,
               out.ind = out.ind) )
}





# ===============================================================
#' @name predic_accuracy
#' @title Return Levels with nonstationarity
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute return levels plot of nonstationary model with the data (in years)
#'
#' @param data numeric vector containing the GEV block-maxima
#' @param gev_nstatio Nonstationary GEV fitted model of class \code{gev.fit}
#' (from package \code{ismev})
#' @param t Maximum time period for which the return levels are considered
#' @param m Return period
#'
#' @return The return levels for the considered time period (t)
#' @examples
#' rl_10_lin <- return.lvl.nstatio(max_years$df$Year,
#'
# DIC Function
#' @rdname predaccur
#' @export
'dic' <- function(out, vals) {
  pm <- colMeans(out) ;   pmv <- log_post1(pm[1], pm[2], pm[3], pm[4], data)
  pmv <- sum(pmv, na.rm = TRUE) ;   vec1 <- rowSums(vals, na.rm = TRUE)
  2*pmv - 4*mean(vec1)
}
#' @rdname predaccur
#' @export
'waic' <- function(vals) {
  vec1 <- log(colMeans(exp(vals))) ;    vec2 <- colMeans(vals)
  sum(2*vec1 - 4*vec2, na.rm = TRUE)
}



# ===============================================================
#' @export crossval.bayes
#' @title Return Levels with nonstationarity
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute return levels plot of nonstationary model with the data (in years)
#'
#' @param data numeric vector containing the GEV block-maxima
#' @param gev_nstatio Nonstationary GEV fitted model of class \code{gev.fit}
#' (from package \code{ismev})
#' @param t Maximum time period for which the return levels are considered
#' @param m Return period
#'
#' @return The return levels for the considered time period (t)
#' @examples
#' rl_10_lin <- return.lvl.nstatio(max_years$df$Year,
#'
"crossval.bayes" <- function(){

}





# ===============================================================
#' @name return_levels_gg
#' @title Return Levels with nonstationarity
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute return levels plot of nonstationary model with the data (in years)
#'
#' @param data numeric vector containing the GEV block-maxima
#' @param gev_nstatio Nonstationary GEV fitted model of class \code{gev.fit}
#' (from package \code{ismev})
#' @param t Maximum time period for which the return levels are considered
#' @param m Return period
#'
#' @return The return levels for the considered time period (t)
#' @examples
#' rl_10_lin <- return.lvl.nstatio(max_years$df$Year,
#'
# POSTERIOR return level plot. Post is the MC generated
# npy is the Number of obs Per Year.
#' @rdname rlfuns
#' @export
"rl.post_gg" <-  function(post, npy, method = c("gev", "gpd"), ci = 0.9, ...) {
  if (method == "gev") npy <- 1

  rps <- c(1/npy + 0.001, 10^(seq(0,4,len=20))[-1])
  p.upper <- 1 - 1/(npy * rps)
  mat <- evdbayes::mc.quant(post = post, p = p.upper, lh = lh)
  mat <- t(apply(mat, 2, quantile, probs = c((1-ci)/2, 0.5, (1+ci)/2)))
  print(mat)
  df <- data.frame('return period' = rps, "TX" = mat  )
  print(df)
  g <- ggplot(df) + geom_line(aes(y = TX.5., x = 'return period'), col= "red") +
    geom_line(aes(y = TX.50., x = 'return period')) +
    geom_line(aes(y = TX.95., x = 'return period'), col = "red") +
    scale_x_log10(breaks = c(1,10,100),labels = c(1,10,100)) +
    theme_piss(...)
  print(g)
  # matplot(rps, mat, log = "x", type = "l",
  #         xlab = xlab, ylab = ylab, lty = lty, col = col, ...)
  return(list(x = rps, y = mat))
}


#' @rdname rlfuns
#' @export
"rl.pred_gg" <- function(post, qlim, npy,
                         method = c("gev", "gpd"), period = 1, ...) {
    if (method == "gev")  npy <- 1

    np <- length(period)
    p.upper <- matrix(0, nrow = 25, ncol = np)
    qnt <- seq(qlim[1], qlim[2], length = 25)

    for(i in 1:25) {
      p <- (qnt[i] - post[,"mu"])/post[,"sigma"]
      p <- ifelse(post[,"xi"],
                  exp( - pmax((1 + post[,"xi"] * p),0)^(-1/post[,"xi"])),
                  exp(-exp(-p)))
      for(j in 1:np)
        p.upper[i,j] <- 1-mean(p^period[j])
    }
    if (method == "gpd")  p <- 1 + log(p)
    if(any(p.upper == 1)) stop("lower q-limit is too small")
    if(any(p.upper == 0)) stop("upper q-limit is too large")

    df <- data.frame(x = 1/(npy * p.upper), y = qnt) ;
    g <- ggplot(df) + geom_line(aes(x = x.1, y = y)) +
      geom_line(aes(x = x.2, y = y)) + geom_line(aes(x = x.3, y = y)) +
      scale_x_log10(breaks = c(1,10,100),labels = c(1,10,100)) +
      theme_piss(...)
   print(g)
    return(list(x = 1/(npy * p.upper ), y = qnt))
}



# ===========================================================================
#' @export pred_post_samples
#' @title Predictive Posterior samples
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#' Compute posterior predictive samples from the obtained model (gibbs.trend)
"pred_post_samples" <- function (from = 1, until = nrow(max_years$df), n_future = 0) {

  tt2 <- ( (max_years$df$Year[from]):(max_years$df$Year[until] + n_future ) -
             mean(max_years$df$Year) )  /  until

  repl2 <- matrix(NA, nrow(gibbs.trend$out.chain), length(tt2))

  for(t in 1:nrow(repl2)) {
    mu <- gibbs.trend$out.chain[t,1] + gibbs.trend$out.chain[t,2] * tt2
    repl2[t,] <- evd::rgev(length(tt2),
                           loc = mu,
                           scale = gibbs.trend$out.chain[t,3],
                           shape = gibbs.trend$out.chain[t,4])
  }
  return(repl2)
}




# ===============================================================
#' @export Pred_Dens_ggPlot
#' @title Predictive Posterior density plots
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}

#' @description
#' Compute a \code{ggplot} of the Density associated with the
#' Posterior Predictive Distribution (PPD)
#' @param year numeric value giving the year for which we want to compute the
#' predictive density
#' @param data_ppd numeric matrix or df of size [n,p] containing the posterior predictive
#'  samples, where \code{n} is the number of simulated samples and \code{p}
#'   is the prediction horizon.
#' @param ylim and \code{xlim} define the grid for the plot.
#'
#' @return a ggplot
#' @examples
#' # For the PPD density plot of the year 2026 :
#' gg4 <- PissoortThesis::Pred_Dens_ggPlot(2026, repl2)
#' # where repl2 is the matrix containing the predictive samples (see Bayes_own_gev.R)
#'
'Pred_Dens_ggPlot' <- function(year, data_ppd ,
                               ylim = c(0.01,.55), xlim = c(27.5,36)) {
  index <- year - 1900 # Retrieve the index from our data series.

  ggplot(data.frame(data_ppd)) +
    stat_density(aes_string(x = paste0("X", index)), geom = "line") +
    labs(x = paste("TX for year", as.character(year))) +
    geom_vline(xintercept = mean(repl2[,index]), col = "blue") +
    geom_vline(xintercept = post.pred2["5%", index], col = "blue") +
    geom_vline(xintercept = post.pred2["95%", index], col = "blue") +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    geom_vline(xintercept = hpd_pred['lower', index], col = "green") +
    geom_vline(xintercept = hpd_pred['upper', index], col = "green") +
    theme_piss()
}



