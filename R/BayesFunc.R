
# ===============================================================
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

# Posterior Density Function
# Compute the log_posterior in a stationary context.
# Be careful to incorporate the fact that the distribution can have finite endpoints.
#' @export
'log_post0' <- function(mu, logsig,xi, data) {
  llhd <- -(gev.nloglik(mu = mu, sig = exp(logsig),
                        xi = xi, data = data))
  lprior <- dnorm(mu, sd = 50, log = TRUE)
  lprior <- lprior + dnorm(logsig, sd = 50, log = TRUE)
  lprior <- lprior + dnorm(xi, sd = 5, log = TRUE)
  lprior + llhd
}



# ===============================================================
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
      geom_line() + theme_piss(18,16, theme_classic()) +
      scale_colour_brewer(name = "chain nr", palette = "Set1") +
      guides(colour = guide_legend(override.aes = list(size= 1.2))),
    ggplot(data, aes(x = iter.chain, y = logsig, col = as.factor(chain.nbr))) +
      geom_line() + theme_piss(18,16, theme_classic()) +
      scale_colour_brewer(name = "chain nr", palette = "Set1") +
      guides(colour = guide_legend(override.aes = list(size= 1.2))),
    ggplot(data, aes(x = iter.chain, y = xi, col = as.factor(chain.nbr))) +
      geom_line() + theme_piss(18,16, theme_classic()) +
      scale_colour_brewer(name = "chain nr", palette = "Set1") +
      guides(colour = guide_legend(override.aes = list(size= 1.2))),
    ncol = 1, top = textGrob(title,
                   gp = gpar(col ="darkolivegreen4",
                             fontsize = 25, font = 4))
  )
}

# ===============================================================
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
# start represents the starting value of the generated chain,
# must explore different ones, typically take the MPLE
# varmat.prop is the variance of the proposal.
#' @export
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
# start is the starting point of the algo. Again, use several ones
# propsd value is as proposed by coles but must tune this value.
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
# Compute Log-posterior  of the model in a nonstationary context.
# By default, we consider the simple linear trend model.
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


# propsd must be weel chosen (e.g. Trial-and-error method)
# start must be a list containing a number of different starting values
#and this number determines the number of chains generated.
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

   print(paste("time is ",round((proc.time() - time)[3], 5), " sec"))
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
  mat <- mc.quant(post = post, p = p.upper, lh = lh)
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
"rl.pred_gg" <- function(post, qlim, npy, method = c("gev", "gpd"), period = 1, ...) {
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
