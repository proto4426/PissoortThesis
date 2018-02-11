'tsDiagGamm' <- function(x, timevar, observed, f = 0.3, type = "normalized") {
  resi <- resid(x$lme, type = type) ;   fits <- fitted(x$lme)
  on.exit(layout(1))
  layout(matrix(1:6, ncol = 3, byrow = TRUE))
  plot(resi ~ fits, ylab = "Normalized Residuals",
       xlab = "Fitted Values", main = "Fitted vs. Residuals")
  lines(lowess(x = fits, y = resi, f = f), col = "blue",
        lwd = 2)
  plot(resi ~ timevar, ylab = "Normalized Residuals",
       xlab = "Time", main = "Time series of residuals")
  lines(lowess(x = timevar, y = resi, f = f), col = "blue", lwd = 2)
  plot(observed ~ fits, ylab = "Observed",
       xlab = "Fitted Values", main = "Fitted vs. Observed",
       type = "n")
  abline(a = 0, b = 1, col = "red")  ;   points(observed ~ fits)
  lines(lowess(x = fits, y = observed, f = f), col = "blue",
        lwd = 2)
  hist(resi, freq = FALSE, xlab = "Normalized Residuals")
  qqnorm(resi) ;    qqline(resi)
  acf(resi, main = "ACF of Residuals")
}

  # Or see this above
# 'Deriv' <- function(mod, n = 200, eps = 1e-7, newdata) {
#   m.terms <- attr(terms(mod), "term.labels")
#   if(missing(newdata)) {
#     newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
#                    function(x) seq(min(x), max(x), length = n))
#     names(newD) <- m.terms
#   } else newD <- newdata
#   X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
#   newD <- newD + eps
#   X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
#   Xp <- (X1 - X0) / eps  ;  Xp.r <- NROW(Xp) ;    Xp.c <- NCOL(Xp)
#   ## dims of bs
#   bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
#   # number of smooth terms
#   t.labs <- attr(mod$terms, "term.labels")  ;   nt <- length(t.labs)
#   ## list to hold the derivatives
#   lD <- vector(mode = "list", length = nt) ;   names(lD) <- t.labs
#   for(i in seq_len(nt)) {
#     Xi <- Xp * 0  ;    want <- grep(t.labs[i], colnames(X1))
#     Xi[, want] <- Xp[, want] ;     df <- Xi %*% coef(mod)
#     df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
#     lD[[i]] <- list(deriv = df, se.deriv = df.sd)
#   }
#   class(lD) <- "Deriv"  ;  lD$gamModel <- mod
#   lD$eps <- eps   ;   lD$eval <- newD - eps
#   return(lD)
# }


'simulate.gamm' <- function(object, nsim = 1, seed = NULL, newdata,
                            freq = FALSE, unconditional = FALSE, ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  if (missing(newdata)) {
    newdata <- object$gam$model
  }

  rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m * n), m, n))
  }

  Rbeta <- rmvn(n = nsim,
                mu = coef(object$gam),
                sig = vcov(object$gam, freq = freq,
                           unconditional = unconditional))
  Xp <- predict(object$gam, newdata = newdata, type = "lpmatrix")
  sims <- Xp %*% t(Rbeta)
  sims
}


'derivSimulCI' <- function(mod, n = 200, eps = 1e-7, newdata, term,
                           samples = 10000) {
  stopifnot(require("MASS"))
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x) - (2*eps), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  newDF <- data.frame(newD) ## needs to be a data frame for predict
  X0 <- predict(mod, newDF, type = "lpmatrix")
  newDF <- newDF + eps
  X1 <- predict(mod, newDF, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  ## sample draws from the posterior distribution of model coefficients
  Rbeta <- t(mvrnorm(n = samples, coef(mod), vcov(mod)))
  ## loop over the terms
  for(i in seq_len(nt)) {
    want <- grep(t.labs[i], colnames(X1))
    lD[[i]] <- list(deriv = Xp[, want] %*% coef(mod)[want],
                    simulations = Xp[, want] %*% Rbeta[want, ])
  }
  class(lD) <- "derivSimulCI"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}



################################################
## Functions for derivatives of GAM(M) models ##
################################################

#  newdata = points along x at which we wish to evaluate the derivative
#  eps  = the distance along x we nudge the points to give us locations p′.
'Deriv' <- function(mod, n = 200, eps = 1e-7, newdata, term) {
  #browser()
  if(inherits(mod, "gamm"))  mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else  newD <- newdata
  newDF <- data.frame(newD) ## needs to be a data frame for predict
  X0 <- predict(mod, newDF, type = "lpmatrix") # lpmatrix forces to return a matrix
  newDF <- newDF + eps
  X1 <- predict(mod, newDF, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)  ;   Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs

  # Loop :Have to work separately with columns of the lpmatrix that relate to each spline
  for(i in seq_len(nt)) {
    Xi <- Xp * 0   # Just to create the matrix of same shape
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5 # SD for the entire spline and not for
       #each of the basis functions that comprise the spline
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod  ;  lD$eps <- eps  ;  lD$eval <- newD - eps
  lD ##return
}

'confint.Deriv' <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else { ## how many attempts to get this right!?!?
    ##term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- df.residual(object$gamModel)
  tVal <- qt(1 - (alpha/2), residual.df)
  ##for(i in term.labs[term]) {
  for(i in term) {
    upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
    lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
    res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}


#plows through the points at which we evaluated the derivative and looks for
#locations where the pointwise 1−α confidence interval doesn’t include zero
'signifD' <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  # Returns list with two components incr and decr which contain the
  #locations where the estimated derivative is positive or negative,
  #respectively, and zero is not contained in the confidence interval.
  list(incr = incr, decr = decr)
}


# S3 method which plots the estimated derivatives and associated confidence intervals.
# displayed as lines or a solid polygon via argument polygon
# sizer = TRUE will colour increasing parts of the spline in blue and in red for the decreasing parts
'plot.Deriv' <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term, eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, main, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term)) term <- term.labs
  else  term <- term.labs[match(term, term.labs)]
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  tVal <- qt(1 - (alpha/2), df.residual(x$gamModel))
  if(missing(ylab))  ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  CI <- confint(x, term = term)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  for(i in term) {
    upr <- CI[[i]]$upper
    lwr <- CI[[i]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,i], x[[i]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,i], rev(x$eval[,i])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,i], upr, lty = "dashed")
      lines(x$eval[,i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}


'plot.derivSimulCI' <- function(x, alpha = 0.05, polygon = TRUE,
                              sizer = FALSE, term,
                              eval = 0, lwd = 3,
                              col = "lightgrey", border = col,
                              ylab, xlab, main, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else {
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  ciFUN <- function(x, alpha) {
    ahalf <- alpha / 2
    apply(x$simulations, 1, quantile, probs = c(ahalf, 1 - ahalf))
  }
  CI <- lapply(x[seq_len(l)], ciFUN, alpha = alpha)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  on.exit(layout(1))
  for(i in term) {
    lwr <- CI[[i]][1,]
    upr <- CI[[i]][2,]
    ylim <- range(upr, lwr)
    plot(x$eval[,i], x[[i]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,i], rev(x$eval[,i])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,i], upr, lty = "dashed")
      lines(x$eval[,i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
    }
  }
  invisible(x)
}





