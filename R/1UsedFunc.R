# ===============================================================
#' @export theme_piss
#' @title Homogeneous theme for ggplots
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#' Theme function for our builded \code{ggplots}.
#' Useful get coherent and similar colours, themes,... for all plots.
#' Also useful to decrease the number of
#'
#' @param size_p Size of the plot's title.
#' @param size_c Size of the axis' title.
#' @param size_l Size of the legend's title.
#' @param theme ggplot's theme for the plot. Set it to NULL if you
#' want put your own theme construction in the "..."
#'
#' @return A personalized ggplot2 theme object to add to every builded plots.
#' @details
#' This function is useful to decrease the amount of code for each ggplots
#' generated as this thesis will use exclusive \code{ggplot2} for the plots.
#' Values of other parameters such as colours could be changed inside the function.
#' @examples
#' # To generate the PP-plot ! See code for more details
#' ggplot(data = data.frame(empirical,model_est), aes(x=empirical,y=model_est)) +
#' geom_point(shape = 1, col = "#33666C") + geom_abline(intercept=0,slope=1,col="red") +
#' theme_piss() +  ggtitle("Probability plot")
#'
"theme_piss" <- function(size_p = 18, size_c = 14,
                         size_l = 12, theme = theme_bw(), ...){

    text <- function(size, ...) element_text(size = size, colour = "#33666C",
                                             face="bold", ...)
    theme +
      theme(plot.title = text(size_p, hjust = 0.5),
            axis.title = text(size_c),
            legend.title = text(size_l)
      ) +
      theme(legend.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = "white"), ...
      )
}




# ===============================================================
#' @export grid_arrange_legend
#' @title Format legend for several ggplots
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#'
#' @return A \code{grid.arrange()} with legend formatted
#' @details
#' This function is useful to decrease the amount of code for each ggplots
#' generated as this thesis will use exclusive \code{ggplot2} for the plots.
#' Values of other parameters such as colours could be changed inside the function.
#' @examples
'grid_arrange_legend' <- function(..., ncol = length(list(...)), nrow = 1,
                                  position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
}


# ===============================================================
#' @export func_season
#' @title Retrieve seasons from months
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#' This function aims to retrieve specific seasons from a month standing
#  in a particular dataset. Helpful to characterize our dataset
#' and the temperatures in terms of seasons
#'
#' @param month of the year from which we want to return the season
#'
#' @return the meteoroglogical season associated to each month
#' @examples
#' # Retrieve seasons based on meteorological seasons
#'
#' TXTN_closed$season <- sapply(TXTN_closed$month, function(x) func_season(x))
#'
'func_season' <- function(month){
  if (month %in% c(12, 01, 02))  return("Winter")
  else if (month %in% c(03, 04, 05))  return( "Spring" )
  else if (month %in% c(06, 07, 08))  return("Summer")
  else  return("Autumn")
}



# ===============================================================
#' @export yearly.extrm
#' @title Yearly extrema retrieval
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#' retrieve the (yearly) extremes for a given list of extrema.
#' The two last arguments specify if we want rather minima (Fun=min, tmp='TN')
#' or maxima. This returns the data AND data.frame with the corresponding year.
#' @param list.y List containing numeric vectors of (daily typically) extrema
#' for each year. Length of \code{list.y} is thus the number of years.
#' @param Fun which function to compute, thus here typically max() or min()
#' @param tmp Specify if we want 'TX' (max. temp.) or 'TN' (min. temp.)
#'
#' @return List of 2 elements : \code{df} which is a data.frame used
#' to compute ggplots. \code{data} is the vector of yearly extrema.
#' @details
#' Easier to handle for ggplots
#' @examples
#' # For maxima (default)
#' max_years <- yearly.extrm()
#' # For minnima
#' min_years <- yearly.extrm(Fun = min, tmp = "TN")
"yearly.extrm" <-
  function(list.y = list_by_years, Fun = max, tmp = 'TX'){

    max_years <- matrix(nrow = 116, ncol = 2)

    for (i in 1:116) {
      max_years[i,1] <- Fun(list.y[[i]][tmp])
    }

    max_years[,2] <- seq(1901, 2016 ,by = 1)
    colnames(max_years) <- c("Max","Year")
    max_data <- max_years[,"Max"]

    if (tmp == "TN") colnames(max_years) <- c("Min", "Year")

    max_frame <- as.data.frame(max_years)
    out <- list(max_data, max_frame)  ;  names(out) <- c("data", "df")
    return(out)
  }


# ===============================================================
#' @export MeanExcess
#' @title Mean Excess function and plots
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#' Computes the mean excess scores \eqn{e_{k,n}} for a numeric vector of data
#' @param data numeric vector of extremes
#' @param plot specify if we want a plot or not
#' @param k specifies if we want the plots as a function of k (=TRUE)
#' or as a function of the order statistics \eqn{X_{n-k,n}} (=FALSE)
#'
#' @return List containing the values of k, the corresponding order statistics
#' \eqn{X_{n-k}} and the mean excess scores associated \code{e_{k,n}}
#' @details
#' Based on the functions provided in \url{http://lstat.kuleuven.be/Wiley/}.
#' If plot=T, the mean excess scores are plotted as a function of k
#' If plot=TRUE and k=TRUE then the mean excess scores are
#' plotted as a function of the order statistics X_n-k,n
#' @examples
#' max_all <- TXTN_closed$TX
#' MeanExcess(max_all, plot = T, k =1000)
#' @references
#' Jan Beirlant, Yuri Goegebeur, Johan Segers, and Jozef Teugels.
#' Statistics of Extremes: Theory and Applications. John Wiley & Sons, March 2006.
#' Section \eqn{\boldsymbol{1.2.2}}
"MeanExcess" <- function(data, plot=FALSE, k=FALSE) {
  X <- sort(data)  ;  n <- length(X)  ;  e <- numeric(n)  ;   K <- 1:(n-1)
  ### estimation of mean excess scores
  for (k in 1:(n-1)) {
    e[k] <- (1/k)*sum(X[n:(n-k+1)]) - X[n-k]
  }
  ### plots if TRUE
  if ( plot ){
    if ( k ) {		### as function of k
      plot(K,e[K], type = "b", ylab = "e_k,n", xlab = "k",
           main = "Mean excess plot")
    }
    else { 	   	### as function of order statistics X_n-k,n
      plot(X[n-K], e[K], type = "b", ylab = "e_k,n",
           xlab = "X_n-k,n", main = "Mean excess plot")
    }
  }
  list(k = K, X = X[n-K], e=e[K])
}

# ===============================================================
#' @export genHill
#' @title Generalized Hill estimates
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#'  Computes the Hill estimates of gamma (Section 5.2.3)
#' for a numeric vector of observations (data) and as a function of k
#' If plot=TRUE then the estimates are plotted as a function of k
#' If add=TRUE then the estimates are added to an existing plot
#' @param data numeric vector of extremes
#' @return list of K and \eqn{\gamma} (Zipf estimates)
#' @details
#' Based on the functions provided in \url{http://lstat.kuleuven.be/Wiley/}.
#' Not used as it holds only for fat-tailed distributions (\eqn{\xi>0})
"genHill" <- function(data, gamma, plot=FALSE, add=FALSE, ...) {
  X <- sort(data)  ;   n <- length(X)  ;   UH.scores <- numeric(n)
  Hill <- numeric(n)     ;     K <- 1:(n-1)
  ### Hill estimates
  for (i in 1:(n-1)) {
    UH.scores[i] <- X[n-i]*gamma[i]
  }
  for (k in 1:(n-2)) {
    Hill[k] <- sum(log(UH.scores[1:k])-log(UH.scores[k+1]))/k
    print(k)
  }
  ### plots if = TRUE
  if (plot || add){
    if ( !add ) {   	### plot estimates
      plot(K, Hill[K], type="l", ylab="gamma", xlab="k",
           main="Estimates of extreme value index", ...)
    }
    else { 			### adds estimates to existing plot
      lines(K, Hill[K], ...)
    }
  } ### output list with values of k and corresponding Zipf estimates
  list(k=K, gamma=Hill[K])
}


# To compute the plots of both ACF and PACF
# ===============================================================
#' @export acfpacf_plot
#' @title (P)ACF plots
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description Produce plots of both ACF and PACF for a
#' given maximum number of lags.
#'
#' @param data numeric ector containing
#' @param laxgmax numeric value containing the maximum number of lags to display
#' @param dotsACF more arguments to pass to the ACF function
#' @param dotsACF more arguments to pass to the ACF function
#' @return Page of both ACF and PACF until lagmax
#' @details Reduce length of code as this operation should be made several times.
#' This will be implemented in \code{ggplot} if shown in the finalo thesis.
#' @examples
#' acfpacf_plot(max_years$data)
'acfpacf_plot' <- function(data, lagmax = length(data), mfrow = c(2, 1),
                           main1 = "ACF", main2 = "PACF"){
  par(mfrow = mfrow)
  acf(data, lag.max = lagmax, main = main1)
  pacf(data, lag.max = lagmax, main = main2)
  # ACF <- acf(data, plot = FALSE)
  # ACF <- setNames(data.frame(unclass(ACF)[c("acf", "lag")]), c("ACF","Lag"))
  # g1 <- ggplot(ACF, aes(x = Lag, y = ACF)) +
  #   geom_hline(aes(yintercept = 0)) +
  #   geom_segment(mapping = aes(xend = Lag, yend = 0))
  #
  #
  # PACF <- pacf(data, plot = FALSE)
  # PACF <- setNames(data.frame(unclass(PACF)[c("pacf", "lag")]), c("PACF","Lag"))
  # g2 <- ggplot(ACF, aes(x = Lag, y = ACF)) +
  #   geom_hline(aes(yintercept = 0)) +
  #   geom_segment(mapping = aes(xend = Lag, yend = 0))
  #
  # grid.arrange(g1,g2, nrow=2)
}


#' @export gev.fit2
#' @title gev.fit() from ismev redefined
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description Redifined function from ismev because there were errors when
#' computing the cubic model for the location parameter in a nonstationary model.
#' We needed to insert a parameter allowing to control the solerance when needed
#' to solve(x$hessian).
#' @param solve.tol and those found in ismev : ?ismev::gev.Fit
#' @examples
#' ## 'Cubic' model ?
#' ti3 <- matrix(ncol = 3, nrow = length(max_years$data))
#' ti3[ ,1] <- seq(1, length(max_years$data), 1)
#' ti3[ ,2] <- (ti3[,1])^2
#' ti3[ ,3] <- (ti3[,1])^3
#' gev_nonstatio3 <- gev.fit(max_years$data/1000, ydat = ti3, mul = c(1, 2, 3))
#' gev_nonstatio3 <- PissoortThesis::gev.fit2(max_years$data, ydat = ti3,
#'                                           mul = c(1, 2, 3),
#'                                           browser = T, solve.tol = 1e-25)
#' # System is singular if we do not scale the data.  But if we scale, the are still
#' #NaNs produced in the covariance matrix.
"gev.fit2"<- function(xdat, ydat = NULL, mul = NULL, sigl = NULL,
                     shl = NULL, mulink = identity, siglink = identity,
                     shlink = identity, muinit = NULL, siginit = NULL,
                     shinit = NULL, show = TRUE, method = "Nelder-Mead",
                     maxit = 10000, browser = F, solve.tol = 1e-18, ...){
  if(browser)  browser()
  z <- list()
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  z$trans <- FALSE	# if maximization fails, could try
  # changing in1 and in2 which are
  # initial values for minimization routine
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
  if(is.null(mul)) {
    mumat <- as.matrix(rep(1, length(xdat)))
    if( is.null( muinit)) muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
    if( is.null( muinit)) muinit <- c(in1, rep(0, length(mul)))
  }
  if(is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdat)))
    if( is.null( siginit)) siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
    if( is.null( siginit)) siginit <- c(in2, rep(0, length(sigl)))
  }
  if(is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if( is.null( shinit)) shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if( is.null( shinit)) shinit <- c(0.1, rep(0, length(shl)))
  }
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  init <- c(muinit, siginit, shinit)
  gev.lik <- function(a) {
    # computes neg log lik of gev model
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    y <- (xdat - mu)/sc
    y <- 1 + xi * y
    if(any(y <= 0) || any(sc <= 0)) return(10^6)
    sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 1))
  }
  x <- optim(init, gev.lik, hessian = TRUE, method = method,
             control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  mu <- mulink(mumat %*% (x$par[1:npmu]))
  sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  z$nllh <- x$value
  z$data <- xdat
  if(z$trans) {
    z$data <-  - log(as.vector((1 + (xi * (xdat - mu))/sc)^(
      -1/xi)))
  }
  z$mle <- x$par
             z$cov <- solve(x$hessian, tol = solve.tol)
  z$se <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, xi)
  if(show) {
    if(z$trans)
      print(z[c(2, 3, 4)])
    else print(z[4])
    if(!z$conv)
      print(z[c(5, 7, 9)])
  }
  class( z) <- "gev.fit"   ;  invisible(z)
}



# ===============================================================
#' @export gev.proflik
#' @title GEV profile neg profile likelihood
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description Computes the negative profile log-likelihood
#' (credits to E. Gilleland and ismev package)
#'
#' @param a is a two dimensionnal vector containing the scale
#'  and the shape parameter respectively
#' @return the negative profile ikelihood
#' @examples
#' for(i in 1:nint) {
#' xp <- x[i]
#' opt <- optim(sol, gev.proflik)
#' sol <- opt$par ; v[i] <- opt$value
#' }
'gev.proflik' <- function(a) {
  if (abs(a[2]) < 10^(-6)) {
    mu <- xp + a[1] * log(-log(1 - p))
    y <- (z$data - mu)/a[1]
    if(is.infinite(mu) || a[1] <= 0)   l <- 10^6
    else l <- length(y) * log(a[1]) + sum(exp(-y)) + sum(y)
  }
  else {
    mu <- xp - a[1]/a[2] * (( - log(1 - p))^( - a[2]) - 1)
    y <- (z$data - mu)/a[1]
    y <- 1 + a[2] * y
    if(is.infinite(mu) || a[1] <= 0 || any(y <= 0))   l <- 10^6
    else   l <- length(y) * log(a[1]) + sum(y^(-1/a[2])) +
      sum(log(y)) * (1/a[2] + 1)
  }
  l
}



# Return Levels with ggplot !
# ===============================================================
#' @name rl_piss
#' @aliases gev.rl.gradient
#' @aliases rl_piss
#' @title Return Levels with ggplot2
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description \code{rl_piss()} is designed to compute plots from \code{ggplot2}
#' package and \code{gev.rl.gradient()} is used inside \code{rl_piss()} to compute
#' the gradient estimates.
#'
#' @param est numeric vector containing the GEV estimates of the 3 parameters.
#' @param p the given probability of exceeding the threshold in a year. Only used
#' through \code{gev.rl.gradient} function.
#' @param mat Covariance matrix of the estimates
#' @param data numeric vector containing the GEV block-maxima
#' @param theme ggplot2's theme used for the plots
#'
#' @return Corresponding plot of the return levels
#' @details
#' See \code{ismev} package details for more in-depth informations
#' @examples
#' rl_piss(gev_tx$mle, gev_tx$cov, gev_tx$data)
#' @references Stuart Coles. An Introduction to Statistical Modeling of
#'  Extreme Values. Springer Series in Statistics. Springer London, London, 2001.
#' @rdname rlfuns
#' @export
"gev.rl.gradient" <- function (est, p) {
  scale <- est[2]  ;   shape <- est[3]
  if (shape < 0)  zero.p <- p == 0
  else   zero.p <- logical(length(p))
  out <- matrix(NA, nrow = 3, ncol = length(p))
  out[1, ] <- 1
  if (any(zero.p)) {
    out[2, zero.p & !is.na(zero.p)] <- rep(-shape^(-1),
                                           sum(zero.p, na.rm = TRUE))
    out[3, zero.p & !is.na(zero.p)] <- rep(scale * (shape^(-2)),
                                           sum(zero.p, na.rm = TRUE))
  }
  if (any(!zero.p)) {
    yp <- -log(1 - p[!zero.p])
    out[2, !zero.p] <- -shape^(-1) * (1 - yp^(-shape))
    out[3, !zero.p] <- scale * (shape^(-2)) * (1 - yp^(-shape)) -
      scale * shape^(-1) * yp^(-shape) * log(yp)
  }
  return(out)
}

#' @rdname rlfuns
#' @export
"rl_piss" <- function(est, mat, data, thm = theme_piss()){
  # eps <- 1e-006   ;   a <- est ;
  # a1 <- a ;   a2 <- a  ;  a3 <- a
  # a1[1] <- a[1] + eps   ;  a2[2] <- a[2] + eps  ;  a3[3] <- a[3] + eps
  f <- c(seq(0.01, 0.999, length = length(data)))
  q <- gevq(est, 1 - f)
  d <- t( gev.rl.gradient( est = est, p=1-f))
  v <- apply(d, 1, ismev::q.form, m = mat)

  df <- data.frame(y = -1/log(f), q = q,
                   upp = q + 1.96 * sqrt(v), low = q - 1.96 * sqrt(v),
                   point = -1/log((1:length(data))/(length(data) + 1)),
                   pdat = sort(data))
  # Plot it
  # g <- ggplot(df) + coord_cartesian(xlim = c(0.1, 1000)) +
  #   geom_line(aes(x = y, y = q), col = "#33666C") +
  #   geom_line(aes( x = y, y = low), col = "red") +
  #   geom_line(aes (x = y, y = upp), col = "red") +
  #   scale_x_log10(breaks = c(0, 1,2, 5, 10,100, 1000),labels = c(0, 1, 2, 5, 10,100, 1000)) +
  #   labs(title = " Return Level plot", x = "Year (log scale)", y = "Quantile") +
  #   geom_point(aes(x = point, y = pdat), col = "#33666C", shape = 1) +
  #   thm
  # g
  return(df = df)
}



# ===============================================================
#' @export return.lvl.nstatio
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
#' @param plot returns a ggplot of the return levels
#' @param out outputs the return levels (TRUE)
#' @param ...  other parmaters to be passed to geom_point() of \code{ggplot}
#' @return  a ggplot in the \code{g} object
#'  and  The return levels for the considered time period (t) in \code{rl}.
#' @examples
#' rl_10_lin <- return.lvl.nstatio(max_years$df$Year,
#'                                 gev_nonstatio, t = 500, m = 10)
'return.lvl.nstatio' <- function( data, gev_nstatio,  start = max(data),
                                  t = 100, m = 10,
                                 plot = T, out = F, ... ){
    y_m <- -(1 / log(1 - 1/m))

    t <- seq(start, start + t, 1)

    rl_m <- (gev_nstatio$mle[1] + gev_nstatio$mle[2] *
               ( t-min(data) )) +
      (gev_nstatio$mle[3] / gev_nstatio$mle[4]) *
      (y_m^gev_nstatio$mle[4] - 1)

    if(plot){
     g <-  ggplot(data.frame(Return.Levels = rl_m, Years = t)) +
        geom_point(aes(x = Years, y = Return.Levels), ...)
    }
          #scale_x_log10(breaks = c(1,10,100),labels = c(1,10,100))

    if(out) return(rl_m)

  return(list(g = g, rl = rl_m) )
}


# ===============================================================
#' @export diplot_fast
#' @title Dispersion Index Plot for threshld selection
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#' Rearranging diplot() function from \code{POT} package
#' because there were speed problems. SEE DETAILS on the \code{POT}'s package
#' @examples
#' data.di <- data.frame(time = seq(1:length(max_all)), obs = max_all)
#' events <- clust(data.di, u = 26, tim.cond = 20/365, clust.max = TRUE)
#' diplot(events[,-3], u.range = c(28,36) ) # too slow
#' diplot_fast(events[,-3], u.range = c(25,35), nt = 1000 )  # Better

'diplot_fast' <- function (data, u.range, nt = max(200, nrow(data))) {
    data <- na.omit(data) ;  date <- data[, "time"]  ;  samp <- data[, "obs"]
    M <- diff(range(date))
    thresh <- seq(u.range[1], u.range[2], length = nt)
    DI <- NULL ; date <- floor(date) ;  time.rec <- length(date)
    for (u in thresh) {
      nb.occ <- NULL
      idx.excess <- samp > u
      lambda <- sum(idx.excess)/M
      for (year in 1:time.rec){
        nb.occ <- c(nb.occ, sum(idx.excess & (date == year)))
        #print(year)
      }
      DI <- c(DI, var(nb.occ)/lambda)
      #print(u)
    }
    conf_sup <- qchisq(1 - (1 - .95)/2, M - 1)/(M - 1)
    conf_inf <- qchisq((1 - .95)/2, M - 1)/(M - 1)
    main <- "Dispersion Index Plot"  ;
    xlab <- "Threshold" ;    ylab <- "Dispersion Index"
    plot(c(thresh, thresh[1]), c(DI, conf_sup), xlab = xlab,
        ylab = ylab, type = "n", main = main)
    rect(0, conf_inf, 2 * u.range[2], conf_sup, col = "lightgrey", border = FALSE)
    lines(thresh, DI)
    return(invisible(list(thresh = thresh, DI = DI)))
}

# ===============================================================
