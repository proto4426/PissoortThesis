
# ===============================================================
#' @name ciBootFunc
#' @author Antoine Pissoort, \email{antoine.pissoort@@student.uclouvain.be}
#' @aliases mle.ksiFun
#' @aliases sdmle.ksiFun
#' @aliases ciBootFunc
#' @aliases ci2Boot.tFunc
#'
#' @title MLE's GPD estimator for \eqn{\xi}
#' and boostrap evaluations
#'
#' @description
#' Function to compute directly our Maximum Likelihood Estimator (MLE) for
#' the shape parameter of the Generalized Pareto Distribution.
#' It isalso possible to determine bootstrapped confidence interals,
#' either by simple or by double bootstrap
#' @param x A numeric vector of data to be fitted.
#' @param B number of bootstrap replicates for simple bootstrap ci
#' @param B1 number of bootstrap replicates for double bootstrap ci (2)
#' for the first stage
#' @param B1 number of bootstrap replicates for double bootstrap ci (2)
#' for the second stage
#' @param level confidence level for the confidence intervals
#' @return MLE computed for the GPD or the standard deviation of this MLE
#' @details
#' This function is useful to decrease the amount of code and compute directly
#' these estimates.
#' @examples
#' X <- above_ured$exceed_red
#' ksi.chap <- gpd_varu_red$mle[2]
#' sd.ksi.chap <- sdmle.ksiFun(X)
#' c500 <- ciBootFunc(B = 500)   # 9 sec.
#' cc <- ci2Boot.tFunc(B1=100,B2=34)
#' @rdname bootfuns
#' @export
'mle.ksiFun' <- function(x) {
  fit0 <- ismev::gpd.fit(xdat = x, threshold = 0, show = F)
  fit <- fit0$mle[2]   ;   return(fit)
}

#' @rdname bootfuns
#' @export
'sdmle.ksiFun' <- function(x) {
  fit0 <- ismev::gpd.fit(xdat = x, threshold = 0, show = F)
  fit <- fit0$se[2]   ;   return(fit)
}
#' @rdname bootfuns
#' @export
'ciBootFunc' <- function(x , B, level = 0.05, print = T)  {
  root.star <- rep(NA,B)     ;   root.star.t  <- rep(NA,B)
  root.starp <- rep(NA,B)    ;    time <- proc.time()
  # Pre-allocation is very important to speed up the loops !!
  for(b in 1:B){
    xstar <- sample(x, replace = T)
    ksi.star <<- mle.ksiFun(xstar)
    root.starp[b] <- (ksi.star - ksi.chap)
    root.star[b] <- sqrt(n) * root.starp[b]
    sigmastar <<- sdmle.ksiFun(xstar)
    root.star.t[b] <-  root.star[b] / sigmastar
  }
  ci.boot <- c(ksi.chap - quantile(root.star, 1 - level/2, na.rm = T)/sqrt(n),
               ksi.chap - quantile(root.star, level/2, na.rm = T)/sqrt(n))
  ci.tboot <- c(-sd.ksi.chap*quantile(root.star.t, 1-level/2, na.rm = T)/sqrt(n),
                -sd.ksi.chap*quantile(root.star.t, level/2, na.rm = T)/sqrt(n)) +
    ksi.chap
  ci.pboot <- c(quantile(root.starp, level/2, na.rm = T),
                quantile(root.starp, 1 - level/2, na.rm = T)) + ksi.chap
  if (print == T) {
    print(paste(c("time (sec) elapsed  for B=", B, (proc.time() - time)[3])))
    print(paste(c("boot L:","U:"),ci.boot))
    print(paste(c("tboot L:","U:"), ci.tboot))
    print(paste(c("pboot L:","U:"), ci.pboot))
  }
  rootOut <- list(ci.boot, ci.pboot, ci.tboot, root.star.t)
  names(rootOut) <- c("bootbasic","bootp","boott", "boot_t")
  return(rootOut)
}

#' @rdname bootfuns
#' @export

'ci2Boot.tFunc' <- function(x , B1, B2, level = 0.05, print = T) {
  ksi.starb1  <- rep(NA,B1)  ; root.star.tb <- ksi.starb1
  ksi.starb2 <- rep(NA,B2)   ;    time <- proc.time()
  for (b1 in 1:B1) {
    xstarb1 <- sample(x, replace = T)
    ksi.starb1[b1] <- mle.ksiFun(xstarb1)
    for (b2 in 1:B2) {
      xstarb2 <- sample(xstarb1, replace = T)
      ksi.starb2[b2] <-  as.numeric(mle.ksiFun(xstarb2))
    }
    varstar2 <- sum(ksi.starb2^2)/B2 - (sum(ksi.starb2)/B2)^2
    sigmastar2 <- sqrt(n) * sqrt(varstar2)
    root.star.tb[b1] <- sqrt(n) * (ksi.starb1[b1] - ksi.chap) / sigmastar2
  }
  sigmastar <- sqrt(n) * (sum(ksi.starb1^2)/B1 - (sum(ksi.starb1)/B1)^2)^{.5}
  ci.2tboot <- c(-sigmastar*quantile(root.star.tb, 1-level/2)/sqrt(n),
                 -sigmastar*quantile(root.star.tb, level/2)/sqrt(n)) +
    ksi.chap
  if (print == T) {
    cat(" B1=", B1, ",B2=", B2, " : time ")
    print(paste(c(proc.time() - time)[3], "sec"))
    cat("---> double_t.boot L:", ci.2tboot[1], "U:", ci.2tboot[2])
  }
  out <- list(root.star.tb, ci.2tboot) ; names(out) <- c("root2t","cit2t")
  return(out)
}

#' @rdname bootfuns
#' @export
## This first function aims at computing the coverage probabilities
#for a vector of different sample sizes contained in n.vec., for simple ci
## M= is the number of MC replications and assumes = B is reasonable
"cov.ci1.vec" <- function (n.vec = n.vec1, M, level = .05) {
  cov <- list(length(n.vec))   ;  covcat <- c()
  ci1b <- matrix(NA,M,6)  # matrix with 2 cols to store length of ci
  t <- proc.time()   ;  timevec <- c()
  for (i in seq_along(n.vec)){
    n <- n.vec[i]    ; nb <- 0
    for (m in 1:M) {
      x.new <- rgpd(n, shape = ksi.chap, scale = sigma, loc = mu)
      ksi.chap <- mle.ksiFun(x.new)
      sd.ksi.chap <- sdmle.ksiFun(x.new)
      # We compute the MC loop and the bootstrap ci replicates  M times both
      ci.eval <<- ciBootFunc(x = x.new, B = M, level = level, print = F)

      x.fevd <- fevd(x.new, type = "GP", threshold = 0)
      ci.prof <- ci.fevd(x.fevd, type = "parameter", which.par = 2,
                         method = "proflik")  #; print(ci.prof)
      nb <- nb + c(indici(ci.eval$bootbasic),
                   indici(ci.eval$boott), indici(ci.eval$bootp),
                   ifelse(ksi >= ci.prof[[1]] & ksi <= ci.prof[[3]], 1, 0))
      # storing for the length of the intervals
      ci1b[m,1] <- ci.eval$bootbasic[2] - ci.eval$bootbasic[1]
      ci1b[m,2] <- ci.eval$bootp[2] - ci.eval$bootp[1]
      ci1b[m,3] <- ci.eval$boott[2] - ci.eval$boott[1]
      ci1b[m,4] <- ci.prof[[3]] - ci.prof[[1]]
    }
    cov[[i]] <- nb/M  ; names(cov[[i]]) <- c("boot","t-boot","p-boot","proflik")
    covcat <- cov[[i]]
    cat("For n=", n.vec[i], "and M=", M, "we have coverages: ","\n")
    print(covcat, max = length(covcat))
    cat("And time", " \n")
    time <- (proc.time() - t)[3]
    print(paste(c(time, "sec")))
    timevec[i] <- time
  }
  # Investigate and compare results with the length of the ci.
  out <- list(cov, timevec, ci1b)
  names(out) <- c("cover", "time", "length")
  return(out)
}

