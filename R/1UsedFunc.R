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
#' @param theme ggplot's theme for the plot.
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
"theme_piss" <-
  function(size_p = 18, size_c = 14, size_l = 12, theme = theme_bw()){
  theme(plot.title = element_text(size = size_p, hjust=0.5,
                                  colour = "#33666C", face="bold"),
        axis.title = element_text(face = "bold", size= size_c,
                                    colour = "#33666C"),
        legend.position = c(.888, .152),
        legend.title = element_text(colour="#33666C", size = size_l, face="bold"),
        legend.background = element_rect(colour = "black"),
        legend.key = element_rect(fill = "white")) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme
}



# This function aims to retrieve specific seasons from month standing
# in a particular dataset
# ===============================================================
#' @export func_season
#' @title Retrieve seasons from months
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @description
#' To characterize our dataset and the temperatures in terms of seasons
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
  if (month %in% c(01, 02, 03))
    return("Winter")
  else if (month %in% c(04, 05, 06))
    return( "Spring" )
  else if (month %in% c(07, 08, 09))
    return("Summer")
  else
    return("Autumn")
}



# ===============================================================
#' @export theme_piss
#' @title Homogeneous theme for ggplots
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
#' @return Page of both ACF and PACF until lagmax
#' @details Reduce length of code as this operation should be made several times.
#' This will be implemented in \code{ggplot} if shown in the finalo thesis.
#' @examples
#' acfpacf_plot(max_years$data)
'acfpacf_plot' <- function(data, lagmax = length(data)){
  par(mfrow = c(2, 1))
  acf(data, lag.max = lagmax)
  pacf(data, lag.max = lagmax)
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
#' @param a numeric vector containing the GEV estimates of the 3 parameters.
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
"gev.rl.gradient" <- function (a, p) {
  scale <- a[2] ;   shape <- a[3]
  if (shape < 0)
    zero.p <- p == 0
  else zero.p <- logical(length(p))
  out <- matrix(NA, nrow = 3, ncol = length(p))
  out[1, ] <- 1
  if (any(zero.p)) {
    out[2, zero.p & !is.na(zero.p)] <- rep(-shape^(-1), sum(zero.p,
                                                            na.rm = TRUE))
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
"rl_piss" <- function(a, mat, data, thm = theme_piss()){
  # eps <- 1e-006
  # a1 <- a ;   a2 <- a  ;  a3 <- a
  # a1[1] <- a[1] + eps   ;  a2[2] <- a[2] + eps  ;  a3[3] <- a[3] + eps
  f <- c(seq(0.01, 0.999, length = length(data)))
  q <- gevq(a, 1 - f)
  d <- t( gev.rl.gradient( a=a, p=1-f))
  v <- apply(d, 1, ismev::q.form, m = mat)

  df <- data.frame(x  = -1/log(f), y = -1/log(f), q = q,
                   upp = q + 1.96 * sqrt(v), low = q - 1.96 * sqrt(v),
                   point = -1/log((1:length(data))/(length(data) + 1)),
                   pdat = sort(data))
  # Plot it
  ggplot(df) + coord_cartesian(xlim = c(0.1, 1000)) +
    geom_line(aes(x = y, y = q), col = "#33666C") +
    geom_line(aes( x = y, y = low), col = "red") +
    geom_line(aes (x = y, y = upp), col = "red") +
    scale_x_log10(breaks = c(1,10,100),labels = c(1,10,100)) +
    labs(title = " Return Level plot", x = "y (log scale)", y = "q") +
    geom_point(aes(x = point, y = pdat), col = "#33666C", shape = 1) +
    thm
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
#'
#' @return The return levels for the considered time period (t)
#' @examples
#' rl_10_lin <- return.lvl.nstatio(max_years$df$Year,
#'                                 gev_nonstatio, t = 500, m = 10)
'return.lvl.nstatio' <-
  function(data, gev_nstatio, t = 100, m = 10 ){
    y_m <- -(1 / log(1 - 1/m))
    t <- seq(max(data), max(data_year) + t, 1)
    rl_m <- (gev_nstatio$mle[1] + gev_nstatio$mle[2] *
               (t-max(max_years$df$Year))) +
      (gev_nstatio$mle[3] / gev_nstatio$mle[4]) *
      (y_m^gev_nstatio$mle[4] - 1)
    g <- ggplot(data.frame(r.lvels = rl_m, years = t)) +
      geom_point(aes(x = years, y = r.lvels)) +
      #scale_x_log10(breaks = c(1,10,100),labels = c(1,10,100)) +
    print(g)
    return(rl_m)
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

'diplot_fast' <-
  function (data, u.range, nt = max(200, nrow(data)))
  {
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
    rect(0, conf_inf, 2 * u.range[2], conf_sup, col = "lightgrey",
         border = FALSE)
    lines(thresh, DI)
    return(invisible(list(thresh = thresh, DI = DI)))
}

# ===============================================================
