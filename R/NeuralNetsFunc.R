# Modifify source function to avoid unesful printings
# ===============================================================
#' @export gevcdn.fit2
#' @title Rebuilding of the function from \code{GEVcdn} package
#'
#' @description
#' this rebuilded version allows us to :
#' 1. Better understand what is going on  and
#' 2. Speed-up function by e.g. reduce unuseful loops or remove printing methods
#' @param x covariate matrix with number of rows equal to the number
#' of samples and number of columns equal to the number of variables.
#' @param y column matrix of target values with number of rows equal7
#' to the # of samples
#' @param sd.norm Weight penalty regularization : sd parameter
#'  for normal distribution prior for the magnitude of input-hiddenlayer
#'   weights; equivalent to weight penalty regularization.#'
#' @return A personalized ggplot2 theme object to add to every builded plots.
#' @details
#' See other function's details in the GEVcdn package
#' @references Cannon, A.J., 2010. A flexible nonlinear modelling framework for nonstationary generalized
#' extreme value analysis in hydroclimatology. Hydrological Processes, 24: 673-685. DOI: 10.1002/hyp.7506
'gevcdn.fit2' <- function (x, y, iter.max = 1000, n.hidden = 2, Th = gevcdn.logistic,
                           fixed = NULL, init.range = c(-0.25, 0.25),
                           scale.min = .Machine$double.eps,
                           beta.p = 3.3, beta.q = 2, sd.norm = Inf, n.trials = 5,
                           method = c("BFGS", "Nelder-Mead"), max.fails = 100, ...) {
  if (!is.matrix(x))  stop("\"x\" must be a matrix")
  if (!is.matrix(y))  stop("\"y\" must be a matrix")
  method <- match.arg(method)
  if (identical(Th, gevcdn.identity))   n.hidden <- 3
  GML <- Inf
  x.min <- apply(x, 2, min) ;   x.max <- apply(x, 2, max)
  x <- sweep(sweep(x, 2, x.min, "-"), 2, x.max - x.min, "/")  # Standardize
  y.min <- min(y) ;   y.max <- max(y)
  y <- (y - y.min)/(y.max - y.min)
  for (i in seq_len(n.trials)) {
    cat(i,"")
    exception <- TRUE  ;    exception.count <- 0
    while (exception) {
      if (identical(names(init.range), c("W1", "W2"))) {
        weights <- unlist(init.range) +
            gevcdn.initialize(x, n.hidden, c(-0.25, 0.25))
      }
      else {
        weights <- gevcdn.initialize(x, n.hidden, init.range)
      }
      fit.cur <- try(suppressWarnings(optim(weights, gevcdn.cost, method = method,
                                            control = list(maxit = iter.max, ...),
                                            x = x, y = y, n.hidden = n.hidden, Th = Th,
                                            fixed = fixed, scale.min = scale.min, beta.p = beta.p,
                                            beta.q = beta.q, sd.norm = sd.norm)), silent = TRUE)
      if (!class(fit.cur) == "try-error")   exception <- fit.cur$value > 1e+308
      if (exception)  exception.count <- exception.count + 1
      if (exception.count == max.fails) {
        stop("exception... check arguments")
      }
    }
    GML.cur <- fit.cur$value
    cat(GML.cur,"")
    if (GML.cur < GML) {
      GML <- GML.cur
      output.cdn <- fit.cur
    }
  }
  cat("\n")
  weights <- output.cdn$par
  cost <- gevcdn.cost(weights, x, y, n.hidden, Th, fixed, scale.min,
                      beta.p, beta.q, sd.norm)
  w <- gevcdn.reshape(x, weights, n.hidden)
  attr(w, "x.min") <- x.min ;     attr(w, "x.max") <- x.max
  attr(w, "y.min") <- y.min ;     attr(w, "y.max") <- y.max
  attr(w, "Th") <- Th  ;          attr(w, "fixed") <- fixed
  attr(w, "scale.min") <- scale.min
  NLL <- attr(cost, "NLL")  ;     penalty <- attr(cost, "penalty")
  attr(w, "GML") <- GML  ;       attr(w, "NLL") <- NLL
  attr(w, "penalty") <- penalty
  if (sd.norm == Inf) {
    if (length(fixed) == 3)  k <- 3
    else {
      if (identical(Th, gevcdn.identity)) {
        k <- (3 - length(fixed)) * (ncol(x) + 1) + length(fixed)
      }
      else  k <- length(weights) - n.hidden * length(fixed)
    }
    n <- nrow(y)
    BIC <- 2 * NLL + k * log(n)  ;    AIC <- 2 * NLL + 2 * k
    AICc <- AIC + (2 * k * (k + 1))/(n - k - 1)
    attr(w, "BIC") <- BIC ;    attr(w, "AIC") <- AIC
    attr(w, "AICc") <- AICc ;    attr(w, "k") <- k
  }
  return(w)
}


### We put functions here, slightly modifying them, just to better understand them

"gevcdn.bag" <- function (x, y, iter.max = 1000, iter.step = 10, n.bootstrap = 30,
            n.hidden = 3, Th = gevcdn.logistic, fixed = NULL,
            init.range = c(-0.25, 0.25), scale.min = .Machine$double.eps,
            beta.p = 3.3, beta.q = 2, sd.norm = Inf,
            method = c("BFGS", "Nelder-Mead"),
            max.fails = 100, silent = TRUE, ...) {

    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")

    method <- match.arg(method)
    if (identical(Th, gevcdn.identity)) n.hidden <- 3
    x.min <- apply(x, 2, min)
    x.max <- apply(x, 2, max)
    x <- sweep(sweep(x, 2, x.min, '-'), 2, x.max - x.min, '/')
    y.min <- min(y)
    y.max <- max(y)
    y <- (y - y.min)/(y.max - y.min)
    w.bootstrap <- vector("list", n.bootstrap)
    for (i in seq_len(n.bootstrap)){
      cat("*** Bootstrap sample", i, "\n")
      cases.in <- sample(nrow(x), replace=TRUE)
      cases.out <- setdiff(1:nrow(x), cases.in)
      x.train <- x[cases.in,,drop=FALSE]
      x.valid <- x[cases.out,,drop=FALSE]
      y.train <- y[cases.in,,drop=FALSE]
      y.valid <- y[cases.out,,drop=FALSE]
      if (identical(names(init.range), c("W1", "W2"))){
        weights <- unlist(init.range) + gevcdn.initialize(x,
                                                          n.hidden, c(-0.25, 0.25))
      } else{
        weights <- gevcdn.initialize(x, n.hidden, init.range)
      }
      cost.best <- cost.train <- cost.valid <- Inf
      iter.best <- iter <- 0
      while(iter < iter.max){
        if (!silent) cat("Iter", iter, ";",
                         sprintf("%.6g", cost.train),
                         ";", sprintf("%.6g", cost.valid), '\n')
        exception <- TRUE
        exception.count <- 0
        while (exception){
          fit <- try(suppressWarnings(optim(weights, gevcdn.cost,
                                            method = method,
                                            control = list(maxit = iter.step, ...),
                                            x = x.train, y = y.train,
                                            n.hidden = n.hidden, Th = Th,
                                            fixed = fixed, scale.min = scale.min,
                                            beta.p = beta.p, beta.q = beta.q,
                                            sd.norm = sd.norm)),
                     silent = TRUE)
          if (!class(fit) == "try-error"){
            exception <- fit$value > 1e+308
          }
          if (exception){
            exception.count <- exception.count + 1
            weights <- gevcdn.initialize(x, n.hidden,
                                         init.range)
            w.best <- gevcdn.reshape(x, weights, n.hidden)
            cost.best <- cost.train <- cost.valid <- Inf
            iter.best <- iter <- 0
          }
          if (exception.count == max.fails){
            stop("exception... check arguments")
          }
        }
        weights <- fit$par
        cost.prev <- cost.train
        cost.train <- gevcdn.cost(weights, x.train, y.train,
                                  n.hidden, Th, fixed, scale.min,
                                  beta.p, beta.q, sd.norm)
        cost.valid <- gevcdn.cost(weights, x.valid, y.valid,
                                  n.hidden, Th, fixed, scale.min,
                                  beta.p, beta.q, sd.norm)
        iter <- iter + iter.step
        if (cost.valid <= cost.best){
          w.best <- gevcdn.reshape(x, weights, n.hidden)
          cost.best <- cost.valid
          iter.best <- iter
        }
        if (abs(cost.train - cost.prev) < .Machine$double.eps){
          cat("local minimum\n")
          break()
        }
      }
      cat("* Best weights at iter", iter.best, ";",
          sprintf("%.6g", cost.best), "\n")
      attr(w.best, "x.min") <- x.min
      attr(w.best, "x.max") <- x.max
      attr(w.best, "y.min") <- y.min
      attr(w.best, "y.max") <- y.max
      attr(w.best, "Th") <- Th
      attr(w.best, "fixed") <- fixed
      attr(w.best, "scale.min") <- scale.min
      attr(w.best, "stopped.training") <- TRUE
      attr(w.best, "cost.valid") <- as.numeric(cost.best)/
        length(cases.out)
      w.bootstrap[[i]] <- w.best
    }
    if (n.bootstrap==1) w.bootstrap <- w.bootstrap[[i]]
    w.bootstrap
}


gevcdn.bootstrap <- function (n.bootstrap, x, y, iter.max = 1000, n.hidden = 2,
            Th = gevcdn.logistic, fixed = NULL,
            init.range = c(-0.25, 0.25), scale.min = .Machine$double.eps,
            beta.p = 3.3, beta.q = 2, sd.norm = Inf, n.trials = 5,
            method = c("BFGS", "Nelder-Mead"),
            boot.method = c("residual", "parametric"),
            init.from.prev = TRUE, max.fails = 100,
            probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
            ...){
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    boot.method <- match.arg(boot.method)
    weights.bootstrap <- parms.bootstrap <- quantiles.bootstrap <- vector("list", n.bootstrap)
    location.bootstrap <- scale.bootstrap <- shape.bootstrap <-
      matrix(NA, ncol=n.bootstrap, nrow=nrow(x))
    for (i in seq_len(n.bootstrap)){
      cat("** Trial", i, "\n")
      if (i==1){
        cat("Fitting model...\n")
        weights <- gevcdn.fit(x, y, iter.max, n.hidden, Th, fixed,
                              init.range, scale.min,
                              beta.p, beta.q, sd.norm,
                              n.trials, method, max.fails, ...)
        parms <- gevcdn.evaluate(x, weights)
        residuals <- (1 + parms[,"shape"]*(y - parms[,"location"])/
                        parms[,"scale"])^(-1/parms[,"shape"])
      }
      if (boot.method=="residual"){
        y.prime <- as.matrix(parms[,"location"] + parms[,"scale"]*
                               (sample(residuals, replace=TRUE)^
                                  (-parms[,"shape"]) - 1)/parms[,"shape"])
      } else if (boot.method=="parametric"){
        y.prime <- y*0
        for(j in seq_len(nrow(y))){
          y.prime[j] <- rgev(1, location = parms[j,"location"],
                             scale = parms[j,"scale"],
                             shape = parms[j,"shape"])
        }
      }
      if (init.from.prev){
        n.trials <- 1
        init.range <- weights
      }
      weights.prime <- gevcdn.fit(x, y.prime, iter.max, n.hidden, Th,
                                  fixed, init.range, scale.min,
                                  beta.p, beta.q, sd.norm,
                                  n.trials, method, max.fails, ...)
      parms.prime <- gevcdn.evaluate(x, weights.prime)
      quantiles.prime <- sapply(probs, qgev,
                                location=parms.prime[,"location"],
                                scale=parms.prime[,"scale"],
                                shape=parms.prime[,"shape"])
      colnames(quantiles.prime) <- probs
      quantiles.bootstrap[[i]] <- quantiles.prime
      weights.bootstrap[[i]] <- weights.prime
      parms.bootstrap[[i]] <- parms.prime
      location.bootstrap[,i] <- parms.prime[,"location"]
      scale.bootstrap[,i] <- parms.prime[,"scale"]
      shape.bootstrap[,i] <- parms.prime[,"shape"]
    }
    list(weights.bootstrap = weights.bootstrap,
         parms.bootstrap = parms.bootstrap,
         location.bootstrap = location.bootstrap,
         scale.bootstrap = scale.bootstrap,
         shape.bootstrap = shape.bootstrap,
         quantiles.bootstrap = quantiles.bootstrap)
}

