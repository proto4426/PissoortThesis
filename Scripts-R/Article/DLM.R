setwd("/home/proto4426/Documents/master_Thesis/PissoortThesis/Scripts-R/Article")

rain_1880 <- read.csv("data/P50_Uccle_1880.csv",
                      sep = " ")
names(rain_1880)[1] <- "Date"

rain_1880$day <- substr(rain_1880$Date,7,8)
rain_1880$month <- as.numeric(substr(rain_1880$Date,5,6))
rain_1880$year <- substr(rain_1880$Date,1,4)

# Retrieve seasons with our function. Based on meteorological seasons
rain_1880$season <- sapply(rain_1880$month,
                           function(x) PissoortThesis::func_season(x))

rain_1901 <- rain_1880[rain_1880$year > 1900,]



### DLM ================================


library(tidyverse)
library(plotly)
library(dlm)

## Random walk + noise
lg <- dlm(m0 = rep(0,2), C0 = 10 * diag(2), FF = matrix(c(1,0),nr=1),
          V = 1.4, GG = matrix(c(1,0,1,1),nr=2), W = diag(c(0,0.2)))

## Time varying DLM
x <- rnorm(100) # covariates
dlr <- dlm(m0 = rep(0,2), C0 = 10 * diag(2), FF = matrix(c(1,0),nr=1),
           V = 1.3, GG = diag(2), W = diag(c(0.4,0.2)),
           JFF = matrix(c(0,1),nr=1), X = as.matrix(x))

unlist(mod)
'ex.book.plot' <- function(ts = rain_1901$RR, dV1 = 15100, dW1 = 0.5 * 1468,
                      dV2 = 15100, dW2 = 0.5 * 1468, brow = F){
  ts <- ts(ts)
  if(brow) browser()
  # mod <- dlmModPoly(order = 1, dV = 15100, dW = 1468)
  # modFilt <- dlmFilter(ts, mod)
  # str(modFilt,1)
  # with(modFilt, dlmSvd2var(U.C[[101]], D.C[101,]))
  mod1 <- dlmModPoly(order = 1, dV = dV1, dW = dW1)
  nileFilt1 <- dlmFilter(ts, mod1)
  plot(window(cbind(ts,nileFilt1$m[-1]),
              start=start(ts)+1), plot.type='s',
       type='o', col=c("grey","green"), lty=c(1,2), xlab="", ylab="Level")
  mod2 <- dlmModPoly(order = 1, dV = dV2, dW = 10* dW2)
  nileFilt2 <- dlmFilter(ts, mod2)
  lines(window(nileFilt2$m,start=start(ts)+1), type='o', col="red", lty=4)
  legend("bottomleft", legend=c("data", "filtered level - model 1",
                                "filtered level - 10x larger signal to noise ratio"),
         col=c("grey", "green", "red"), lty=c(1,2,4), pch=1, bty='n')
  
  return(list(mod1 = nileFilt1, mod2 = nileFilt2))
}
#ex.book.plot()

list_rain_by_years <- split(rain_1901, rain_1901$year)
max_years <- PissoortThesis::yearly.extrm(list.y = list_rain_by_years, tmp = "RR")

nileFilt1 <- ex.book.plot(max_years$df$Max,)$mod1
nileFilt2 <- ex.book.plot(max_years$df$Max,)$mod2

maxy <- ex.book.plot(max_years$df$Max, dV1 = 10000, dW1 = 50000)


'ex.book.ggplot' <- function(ts = rain_1901$RR, dV = 15100, dW = 0.5 * 1468,
                             dV2 = 15100, dW2 = 0.5 * 1468, brow = F){
  
  ts <- ts(ts)
  mod1 <- dlmModPoly(order = 1, dV = dV, dW = dW)
  nileFilt1 <- dlmFilter(ts, mod1)

  mod2 <- dlmModPoly(order = 1, dV = dV2, dW = dW2)
  nileFilt2 <- dlmFilter(ts, mod2)
  
  df1 <- as.data.frame(window(cbind(series = ts, 
                                    filter1 = nileFilt1$m[-1],
                                    filter2 = nileFilt2$m[-1],
                                    obs = 1:length(ts)),
                              start=start(ts)+1))

  cat("var of the filtering distribution of the last-state vector",
      with(nileFilt1, dlmSvd2var(U.C[[length(ts)+1]], D.C[length(ts)+1, ])) ) 
  
  g <- ggplot(df1, aes(x = obs)) + 
    geom_point(aes(y = series)) + geom_line(aes(y = series)) +
    geom_line(aes(y = filter2), col = "red")  + 
    geom_line(aes(y = filter1), col = "green3") + 
    annotate(geom = "text", label = paste("dV2 = ", round(dV2,1), "(dW2=", round(dW2),")"),
             x = 20, y = 65, col = "red" , size = 3) +
    annotate(geom = "text", label = paste("dW = ", round(dW, 1 ), "(dV=", round(dV),")"),
             x = 20, y = 60, col = "green3" , size = 3) +
    PissoortThesis::theme_piss()
  #print(g)
  
  return(list(mod = nileFilt1, g = g))
}

maxy <- ex.book.ggplot(max_years$df$Max, dV = 10000, dW = 5000)
ggplotly(maxy$g + labs(title = "xxxx"))


for (i in 1:12){
  g <- ex.book.ggplot(max_years$df$Max, 
                      dV2 = 1500 * (i^1.5), 
                      dW = 200 * (i^1.5) )$g
  assign(paste0("g", i), g)
}
library(gridExtra)
library(grid)

glist <- as.list(paste0("g", seq(1,12)[1]))
n <- length(glist)  ;   nCol <- floor(sqrt(n))
do.call("grid.arrange", c(lapply(glist, [1]), ncol=nCol))

grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8, g9, g10, g11, g12, ncol = 2)

cowplot::plot_grid(plotlist = glist, ncol = 2)


###### Smoothing 

ts <- ts(max_years$df$Max)
mod1 <- dlmModPoly(order = 1, dV = 15100, dW = 0.5 * 1468)
modFilt <- dlmFilter(ts, mod1)
modSmooth <- dlmSmooth(modFilt)
str(modSmooth,1)
with(modSmooth, drop(dlmSvd2var(U.S[[length(ts)]], D.S[length(ts),])))
with(modFilt, drop(dlmSvd2var(U.C[[length(ts)/2]], D.C[length(ts)/2,])))
with(modSmooth, drop(dlmSvd2var(U.S[[length(ts)/2]], D.S[length(ts)/2,])))

hwid <- qnorm((1-0.95) / 2) *
  sqrt(with(modSmooth, unlist(dlmSvd2var(U.S, D.S))))
smooth <- cbind(modSmooth$s, as.vector(modSmooth$s) + hwid %o% c(1,-1))
plot(cbind(ts, window(smooth, start=start(ts))), plot.type="s",
     col=c("grey", "magenta", "cyan", "cyan"), lty=c(1, 2, 3, 3), type="o",
     ylab="Level", xlab="")
legend("bottomleft", legend=c("data", "smoothed level", "95% probability limits"),
       col=c("grey", "magenta", "cyan"), lty=1:3, pch=1, bty="n")


mod <- dlm(m0 = rep(0,4), C0 = 1e-8 * diag(4),
             FF = matrix(c(1,1,0,0), nr=1),
             V = 1e-3,
             GG = bdiag(matrix(1),
                          matrix(c(-1,-1,-1,1,0,0,0,1,0),nr=3,byrow=TRUE)),
             W = diag(c(771.35, 86.48, 0, 0), nr=4))
modFilt <- dlmFilter(ts, mod)

modSmooth <- dlmSmooth(modFilt)
plot(ts, type="o", xlim=c(1957, 1970), ylim=c(210, 725), ylab="Expenditure")
lines(modFilt$m[,1], col="green", type="o")
lines(modSmooth$s[,1], col="red", type="o")
plot(modSmooth$s[,2], type="o", ylab="Expenditure - Seasonal component")
abline(h=0)


## 2.5. Forecasting 

plot(window(cbind(ts, nileFilt1$f, nileFilt2$f), start=1, end=116+10),
       plot.type="s", type="o", col=c("grey", "green", "red"), lty=c(1,2,4),
       xlab="", ylab="Level")
legend("bottomleft", legend=c("data", paste("one-step-ahead forecast - model", 1:2)),
           col=c("grey", "green", "red"), lty=c(1,2,4), pch=1, bty="n")


fore <- dlmForecast(modFilt, nAhead = 8, sampleNew = 10)
invisible(lapply(fore$newObs, function(x) lines(x, col="grey", lty=2)))
lines(fore$f, col="magenta", type="o")



## 3.2.1  Trend models

# standardized innovations and perform Shapiro-Wilk and Ljung-Box tests.
modLSup <- dlmModPoly(1, dV = 9.465, dW = 0.121)
lSupFilt <- dlmFilter(ts, modLSup)
res <- residuals(lSupFilt, sd=FALSE)
stats::shapiro.test(res)
Box.test(res, lag=20, type="Ljung")


### 4.6 Gibbs Sampling ....

dlmBSample
arms
gdpGibbs
dlmGibbsDIG
dlmGibbsDIGt



#### DLM PACKAGE VIGNETTES #############################

## 1.3 Combining Models 

myMod <- dlmModPoly() + dlmModSeas(4) # stochastic linear trend + 4terly seasonal component

# Outer sum : dimension of each model may be different.

dlmModPoly(dV = 0.2, dW = c(0, 0.5)) %+%
  (dlmModSeas(4, dV = 0, dW = c(0, 0, 0.35)) +
     dlmModPoly(1, dV = 0.1, dW = 0.03))


## 1.4 Time-varying Models 

u <- rnorm(25)
myMod <- dlmModReg(u, dV = 14.5)
myMod$JFF
head(myMod$X)


### 2. MLE 

dlmMLE   ; dlmLL

#random walk plus noise, with unknown system and observation variances. 
#Parametrizing variances on a log scale, to ensure positivity,
buildFun <- function(x)   dlmModPoly(1, dV = exp(x[1]), dW = exp(x[2]))

fit <- dlmMLE(ts, parm = c(0,0), build = buildFun) ;  fit$conv
dlmNile <- buildFun(fit$par)
V(dlmNile)
W(dlmNile)
StructTS(ts)


#take into account a jump #in the flow of the river following the construction 
#f Ashwan dam in 1898. this can be done by
'buildFun' <- function(x) {
  m <- dlmModPoly(1, dV = exp(x[1]))
  m$JW <- matrix(1)
  m$X <- matrix(exp(x[2]), nc = 1, nr = length(Nile))
  j <- which(time(Nile) == 1899)
  m$X[j,1] <- m$X[j,1] * (1 + exp(x[3]))
  return(m)
}
fit <- dlmMLE(ts, parm = c(0,0,0), build = buildFun)  ;  fit$conv

dlmNileJump <- buildFun(fit$par)
V(dlmNileJump)
dlmNileJump$X[c(1, which(time(Nile) == 1899)), 1]

### 3.1 Filtering 

nileJumpFilt <- dlmFilter(Nile, dlmNileJump)
plot(Nile, type = 'o', col = "seagreen")
lines(dropFirst(nileJumpFilt$m), type = "o", pch = 20, col = "brown")

attach(nileJumpFilt)
v <- unlist(dlmSvd2var(U.C, D.C))
pl <- dropFirst(m) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(m) + qnorm(0.95, sd = sqrt(v[-1]))

detach()
lines(pl, lty = 2, col = "brown")
lines(pu, lty = 2, col = "brown")


### 3.2 Smoothing

nileJumpSmooth <- dlmSmooth(nileJumpFilt)
plot(Nile, type = 'o', col = "seagreen")

attach(nileJumpSmooth)
lines(dropFirst(s), type = 'o', pch = 20, col = "brown")
v <- unlist(dlmSvd2var(U.S, D.S))
pl <- dropFirst(s) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(s) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown")
lines(pu, lty = 2, col = "brown")


# lGas example : quarterly seasonal component and a local linear trend,
#in the form of an integrated random walk.

lGas <- log(UKgas)
dlmGas <- dlmModPoly() + dlmModSeas(4)
buildFun <- function(x) {
  diag(W(dlmGas))[2:3] <- exp(x[1:2])
  V(dlmGas) <- exp(x[3])
  return(dlmGas)
}
(fit <- dlmMLE(lGas, parm = rep(0, 3), build = buildFun))$conv

dlmGas <- buildFun(fit$par)
drop(V(dlmGas))
diag(W(dlmGas))[2:3]

gasSmooth <- dlmSmooth(lGas, mod = dlmGas)
x <- cbind(lGas, dropFirst(gasSmooth$s[,c(1,3)]))
colnames(x) <- c("Gas", "Trend", "Seasonal")
plot(x, type = 'o', main = "UK Gas Consumption")


### 3.3. Forecasting

# Means and variances of future states and observations are returned
#in a list as components a, R, f, and Q.
gasFilt <- dlmFilter(lGas, mod = dlmGas)
gasFore <- dlmForecast(gasFilt, nAhead = 20)
sqrtR <- sapply(gasFore$R, function(x) sqrt(x[1,1]))
pl <- gasFore$a[,1] + qnorm(0.05, sd = sqrtR)
pu <- gasFore$a[,1] + qnorm(0.95, sd = sqrtR)
x <- ts.union(window(lGas, start = c(1982, 1)),
              window(gasSmooth$s[,1], start = c(1982, 1)),
              gasFore$a[,1], pl, pu)
plot(x, plot.type = "single", type = 'o', pch = c(1, 0, 20, 3, 3),
     col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"),
     ylab = "Log gas consumption")
legend("bottomright", legend = c("Observed", "Smoothed (deseasonalized)",
                                 "Forecasted level", "90% probability limit"),
       bty = 'n', pch = c(1, 0, 20, 3, 3), lty = 1,
       col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"))


#### 4. Bayesian analysis of DLM

## 4.1. Forward filtering backward sampling
plot(Nile, type = 'o', col = "seagreen")
nileFilt <- dlmFilter(Nile, dlmNile)
for (i in 1:10) # 10 simulated "true" levels
    lines(dropFirst(dlmBSample(nileFilt)), col = "brown")

## 4.2. Adaptive rejection Metropolis sampling
lmixnorm <- function(x, weights, means, sds) {
    log(crossprod(weights, exp(-0.5 * ((x - means) / sds)^2
                                 - log(sds))))
}
y <- arms(0, myldens = lmixnorm, indFunc = function(x,...) 
             (x > (-100)) * (x < 100), n = 5000, weights = c(1, 3, 2),
            means = c(-10, 0, 10), sds = c(7, 5, 2))
summary(y)

library(MASS)
truehist(y, prob = TRUE, ylim = c(0, 0.08), bty = 'o')
curve(colSums(c(1, 3, 2) / 6 * dnorm(matrix(x, 3, length(x), TRUE),
                      mean = c(-10, 0, 10), sd = c(7, 5, 2))), add = T)
legend(-25, 0.07, "True density", lty = 1, bty = 'n')


## ARMS Examples ========================================================
#### ==> Warning: running the examples may take a few minutes! <== ####    
## Not run: 
set.seed(4521222)
### Univariate densities
## Unif(-r,r) 
y <- arms(runif(1,-1,1), function(x,r) 1, function(x,r) (x>-r)*(x<r), 5000, r=2)
summary(y); hist(y,prob=TRUE,main="Unif(-r,r); r=2")
## Normal(mean,1)
norldens <- function(x,mean) -(x-mean)^2/2 
y <- arms(runif(1,3,17), norldens, function(x,mean) ((x-mean)>-7)*((x-mean)<7),
          5000, mean=10)
summary(y); hist(y,prob=TRUE,main="Gaussian(m,1); m=10")
curve(dnorm(x,mean=10),3,17,add=TRUE)
## Exponential(1)
y <- arms(5, function(x) -x, function(x) (x>0)*(x<70), 5000)
summary(y); hist(y,prob=TRUE,main="Exponential(1)")
curve(exp(-x),0,8,add=TRUE)
## Gamma(4.5,1) 
y <- arms(runif(1,1e-4,20), function(x) 3.5*log(x)-x,
          function(x) (x>1e-4)*(x<20), 5000)
summary(y); hist(y,prob=TRUE,main="Gamma(4.5,1)")
curve(dgamma(x,shape=4.5,scale=1),1e-4,20,add=TRUE)
## Gamma(0.5,1) (this one is not log-concave)
y <- arms(runif(1,1e-8,10), function(x) -0.5*log(x)-x,
          function(x) (x>1e-8)*(x<10), 5000)
summary(y); hist(y,prob=TRUE,main="Gamma(0.5,1)")
curve(dgamma(x,shape=0.5,scale=1),1e-8,10,add=TRUE)
## Beta(.2,.2) (this one neither)
y <- arms(runif(1), function(x) (0.2-1)*log(x)+(0.2-1)*log(1-x),
          function(x) (x>1e-5)*(x<1-1e-5), 5000)
summary(y); hist(y,prob=TRUE,main="Beta(0.2,0.2)")
curve(dbeta(x,0.2,0.2),1e-5,1-1e-5,add=TRUE)
## Triangular
y <- arms(runif(1,-1,1), function(x) log(1-abs(x)), function(x) abs(x)<1, 5000)     
summary(y); hist(y,prob=TRUE,ylim=c(0,1),main="Triangular")
curve(1-abs(x),-1,1,add=TRUE)
## Multimodal examples (Mixture of normals)
lmixnorm <- function(x,weights,means,sds) {
  log(crossprod(weights, exp(-0.5*((x-means)/sds)^2 - log(sds))))
}
y <- arms(0, lmixnorm, function(x,...) (x>(-100))*(x<100), 5000, weights=c(1,3,2),
          means=c(-10,0,10), sds=c(1.5,3,1.5))
summary(y); hist(y,prob=TRUE,main="Mixture of Normals")
curve(colSums(c(1,3,2)/6*dnorm(matrix(x,3,length(x),byrow=TRUE),c(-10,0,10),c(1.5,3,1.5))),
      par("usr")[1], par("usr")[2], add=TRUE)

### Bivariate densities 
## Bivariate standard normal
y <- arms(c(0,2), function(x) -crossprod(x)/2,
          function(x) (min(x)>-5)*(max(x)<5), 500)
plot(y, main="Bivariate standard normal", asp=1)
## Uniform in the unit square
y <- arms(c(0.2,.6), function(x) 1, function(x) (min(x)>0)*(max(x)<1), 500)
plot(y, main="Uniform in the unit square", asp=1)
polygon(c(0,1,1,0),c(0,0,1,1))
## Uniform in the circle of radius r
y <- arms(c(0.2,0), function(x,...) 1,
          function(x,r2) sum(x^2)<r2, 500, r2=2^2)
plot(y, main="Uniform in the circle of radius r; r=2", asp=1)
curve(-sqrt(4-x^2), -2, 2, add=TRUE)
curve(sqrt(4-x^2), -2, 2, add=TRUE)
## Uniform on the simplex
simp <- function(x) if ( any(x<0) || (sum(x)>1) ) 0 else 1
y <- arms(c(0.2,0.2), function(x) 1, simp, 500)
plot(y, xlim=c(0,1), ylim=c(0,1), main="Uniform in the simplex", asp=1)
polygon(c(0,1,0), c(0,0,1))
## A bimodal distribution (mixture of normals)
bimodal <- function(x) { log(prod(dnorm(x,mean=3))+prod(dnorm(x,mean=-3))) }
y <- arms(c(-2,2), bimodal, function(x) all(x>(-10))*all(x<(10)), 500)
plot(y, main="Mixture of bivariate Normals", asp=1)

## A bivariate distribution with non-convex support
support <- function(x) {
  return( as.numeric( -1 < x[2] && x[2] < 1 && -2 < x[1] &&
                       ( x[1] < 1 || crossprod(x-c(1,0)) < 1 ) ) )
}
Min.log <- log(.Machine$double.xmin) + 10
logf <- function(x) {
  if ( x[1] < 0 ) return(log(1/4))
  else
    if (crossprod(x-c(1,0)) < 1 ) return(log(1/pi))
  return(Min.log)
}
x <- as.matrix(expand.grid(seq(-2.2,2.2,length=40),seq(-1.1,1.1,length=40)))
y <- sapply(1:nrow(x), function(i) support(x[i,]))
plot(x,type='n',asp=1)
points(x[y==1,],pch=1,cex=1,col='green')
z <- arms(c(0,0), logf, support, 1000)
points(z,pch=20,cex=0.5,col='blue')
polygon(c(-2,0,0,-2),c(-1,-1,1,1))
curve(-sqrt(1-(x-1)^2),0,2,add=TRUE)
curve(sqrt(1-(x-1)^2),0,2,add=TRUE)
sum( z[,1] < 0 ) # sampled points in the square
sum( apply(t(z)-c(1,0),2,crossprod) < 1 ) # sampled points in the circle
## End(Not run)
# ===================================


### 4.3. Gibbs sampling : example
dlmGibbsDIG <- function (y, mod, a.y, b.y, a.theta, b.theta, shape.y, rate.y, 
          shape.theta, rate.theta, n.sample = 1, thin = 0, ind, save.states = TRUE, 
          progressBar = interactive()) 
{
  msg1 <- "Either \\"a.y\\" and \\"b.y\\" or \\"shape.y\\" and \\"rate.y\\" must be specified"
  msg2 <- "Unexpected length of \\"shape.y\\" and/or \\"rate.y\\""
  msg3 <- "Unexpected length of \\"a.y\\" and/or \\"b.y\\""
  msg4 <- paste("Either \\"a.theta\\" and \\"b.theta\\" or \\"shape.theta\\"", 
                "and \\"rate.theta\\" must be specified")
  msg5 <- "Unexpected length of \\"shape.theta\\" and/or \\"rate.theta\\""
  msg6 <- "Unexpected length of \\"a.theta\\" and/or \\"b.theta\\""
  msg7 <- "\\"thin\\" must be a nonnegative integer"
  msg8 <- "multivariate observations are not allowed"
  msg9 <- "inadmissible value of \\"ind\\""
  if (NCOL(y) > 1) 
    stop(msg8)
  r <- ncol(mod$FF)
  if (hasArg(ind)) {
    ind <- unique(as.integer(ind))
    s <- 1:r
    if (!all(ind %in% s)) 
      stop(msg9)
    perm <- s[c(ind, s[!(s %in% ind)])]
    FF(mod) <- mod$FF[, perm, drop = FALSE]
    GG(mod) <- mod$GG[perm, perm, drop = FALSE]
    W(mod) <- mod$W[perm, perm, drop = FALSE]
    p <- length(ind)
  }
  else {
    perm <- ind <- 1:r
    p <- r
  }
  nobs <- NROW(y)
  if (is.numeric(thin) && (thin <- as.integer(thin)) >= 0) {
    every <- thin + 1
    mcmc <- n.sample * every
  }
  else stop(msg7)
  if (!hasArg(a.y)) 
    if (!hasArg(shape.y)) 
      stop(msg1)
  else if (!hasArg(rate.y)) 
    stop(msg1)
  else {
    if (!all(c(length(shape.y), length(rate.y)) == 1)) 
      warning(msg2)
  }
  else if (!hasArg(b.y)) 
    stop(msg1)
  else {
    if (!all(c(length(a.y), length(b.y)) == 1)) 
      warning(msg3)
    shape.y <- a.y^2/b.y
    rate.y <- a.y/b.y
  }
  if (!hasArg(a.theta)) 
    if (!hasArg(shape.theta)) 
      stop(msg4)
  else if (!hasArg(rate.theta)) 
    stop(msg4)
  else {
    if (!all(c(length(shape.theta), length(rate.theta)) %in% 
             c(1, p))) 
      warning(msg5)
  }
  else if (!hasArg(b.theta)) 
    stop(msg4)
  else {
    if (!all(c(length(a.theta), length(b.theta)) %in% c(1, 
                                                        p))) 
      warning(msg6)
    shape.theta <- a.theta^2/b.theta
    rate.theta <- a.theta/b.theta
  }
  shape.y <- shape.y + 0.5 * nobs
  shape.theta <- shape.theta + 0.5 * nobs
  shape.theta <- rep(shape.theta, length.out = p)
  rate.theta <- rep(rate.theta, length.out = p)
  theta <- matrix(0, nobs + 1, r)
  if (save.states) 
    gibbsTheta <- array(0, dim = c(nobs + 1, r, n.sample))
  gibbsV <- vector("numeric", n.sample)
  gibbsW <- matrix(0, nrow = n.sample, ncol = p)
  it.save <- 0
  if (progressBar) 
    pb <- txtProgressBar(0, mcmc, style = 3)
  for (it in 1:mcmc) {
    if (progressBar) 
      setTxtProgressBar(pb, it)
    modFilt <- dlmFilter(y, mod, simplify = TRUE)
    theta[] <- dlmBSample(modFilt)
    y.center <- y - tcrossprod(theta[-1, , drop = FALSE], 
                               mod$FF)
    SSy <- drop(crossprod(y.center))
    mod$V[] <- 1/rgamma(1, shape = shape.y, rate = rate.y + 
                          0.5 * SSy)
    theta.center <- theta[-1, , drop = FALSE] - tcrossprod(theta[-(nobs + 
                                                                     1), , drop = FALSE], mod$GG)
    SStheta <- drop(sapply(1:p, function(i) crossprod(theta.center[, 
                                                                   i])))
    SStheta <- colSums((theta[-1, 1:p, drop = FALSE] - tcrossprod(theta[-(nobs + 
                                                                            1), , drop = FALSE], mod$GG)[, 1:p])^2)
    diag(mod$W)[1:p] <- 1/rgamma(p, shape = shape.theta, 
                                 rate = rate.theta + 0.5 * SStheta)
    if (!(it%%every)) {
      it.save <- it.save + 1
      if (save.states) 
        gibbsTheta[, , it.save] <- theta
      gibbsV[it.save] <- diag(mod$V)
      gibbsW[it.save, ] <- diag(mod$W)[1:p]
    }
  }
  colnames(gibbsW) <- paste("W", ind, sep = ".")
  if (progressBar) 
    close(pb)
  if (save.states) 
    return(list(dV = gibbsV, dW = gibbsW, 
                theta = gibbsTheta[, order(perm), , drop = FALSE]))
  else return(list(dV = gibbsV, dW = gibbsW))
}

library(dlm)
library(PissoortThesis)
data(max_years)

outGibbs <- dlmGibbsDIG(max_years$df$Max, dlmModPoly(2) + dlmModSeas(4),
                          a.y = 1, b.y = 1000, a.theta = 1,
                          b.theta = 1000, n.sample = 1100, ind = c(2, 3),
                          save.states = FALSE)
burn <- 100
attach(outGibbs)
dV <- dV[-(1:burn)]
dW <- dW[-(1:burn), ]
detach()
par(mfrow=c(2,3), mar=c(3.1,2.1,2.1,1.1))
plot(dV, type = 'l', xlab = "", ylab = "",
     main = expression(sigma^2))
plot(dW[ , 1], type = 'l', xlab = "", ylab = "",
     main = expression(sigma[beta]^2))
plot(dW[ , 2], type = 'l', xlab = "", ylab = "",
     main = expression(sigma[s]^2))
use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0, use), n = 3); at <- at[at >= from]
plot(ergMean(dV, from), type = 'l', xaxt = 'n',xlab = "", ylab = "")
axis(1, at = at - from, labels = format(at))
plot(ergMean(dW[ , 1], from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at = at - from, labels = format(at))
plot(ergMean(dW[ , 2], from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at = at - from, labels = format(at))

#===================================================================================
library(nimble)

## Building a simple linear state-space model.
## x is latent space, y is observed data
 <- nimbleCode({
  x[1] ~ dnorm(mu_0, 1)
  y[1] ~ dnorm(x[1], 1)
  for(i in 2:t){
    x[i] ~ dnorm(x[i-1] * a + b, 1)
    y[i] ~ dnorm(x[i] * c, 1)
  }
  a ~ dunif(0, 1)
  b ~ dnorm(0, 1)
  c ~ dnorm(1,1)
  mu_0 ~ dnorm(0, 1)
})
## simulate some data
t <- 25; mu_0 <- 1
x <- rnorm(1 ,mu_0, 1)
y <- rnorm(1, x, 1)
a <- 0.5; b <- 1; c <- 1
for(i in 2:t){
  x[i] <- rnorm(1, x[i-1] * a + b, 1)
  y[i] <- rnorm(1, x[i] * c, 1)
}
## build the model
rTimeModel <- nimbleModel(timeModelCode, constants = list(t = t),
                          data <- list(y = y), check = FALSE )
## Set parameter values and compile the model
rTimeModel$a <- 0.5
rTimeModel$b <- 1
rTimeModel$c <- 1
rTimeModel$mu_0 <- 1
cTimeModel <- compileNimble(rTimeModel)


## ============================== http://sbfnk.github.io/mfiidd/pmcmc.html =======


'my_particleFilter' <- function(fitmodel, theta, init.state,
                                data, n.particles, brow=F) {
  # - fitmodel: a fitmodel object
  # - theta: named numeric vector. Values of the parameters for
  #           which the marginal log-likelihood is desired.
  # - init.state: named numeric vector. Initial values of the state variables.
  # - data: data frame. Observation times and observed data.
  # The function returns the value of the marginal log-likelihood
  
  ## Initialisation of the algorithm
  if(brow) browser()
  # Marginal log-likelihood is set to 0 and will be updated 
  #during the filtering steps
  margLogLike <- 0
  # Particle states can be stored in a list
  state.particles  <- rep(list(init.state), n.particles)
  # Weight: initially equal for all the particles 
  # particle weight can be stored in a vector
  weight.particles <- rep(1/n.particles, length = n.particles)
  # Initialise time variable
  current.time <- 0
  
  ## Loop over observation times: resample, propagate, weight
  for(i in seq_len(nrow(data))){
    # Extract next data point (must be a vector)
    data.point <- unlist(data[i, ])
    next.time <- data.point["time"]
    
    # Resample particles according to their weights. 
    # You can use the `sample` function of R
    # (normalisation of the weights is done in the function)
    index.resampled <- sample(x = n.particles,
                              size = n.particles,
                              replace = TRUE,
                              prob = weight.particles)
    state.particles <- state.particles[index.resampled]
    
    ## Loop over particles: propagate and weight
    for(p in 1:n.particles){
      # Extract current state of the particle 
      current.state.particle <- state.particles[[p]]
      
      # Propagate the particle from current observation time 
      # to the next one using the function `fitmodel$simulate`
      traj <- fitmodel$simulate(theta = theta,
                                init.state = current.state.particle,
                                times = c(current.time,next.time))
      
      # Extract state of the model at next observation time
      # Also make sure that model.point is a vector
      model.point <- unlist(traj[2,fitmodel$state.names])
      
      # Weight the particle with the likelihood of the observed 
      # data point using the function `fitmodel$dPointObs`
      weight.particles[p] <-
        fitmodel$dPointObs(data.point = data.point,
                           model.point = model.point,
                           theta = theta)
      # Update state of the p particle
      state.particles[[p]] <- model.point
    }
    # Increment time
    current.time <- next.time
    
    ## Increment the marginal log-likelihood
    # Add the log of the mean of the particles weights
    margLogLike <- margLogLike + log(mean(weight.particles))
  }
  
  ## Return marginal log-likelihood
  return(margLogLike)
}

library(fitR)
data(SEIT4L_stoch)
data(FluTdC1971)

# theta close to the mean posterior estimate of the deterministic SEIT4L
# model
theta <- c(R0 = 7, D_lat = 1, D_inf = 4,
           alpha = 0.5, D_imm = 10, rho = 0.65)
# init state 
init.state <- c(S = 279, E = 0, I = 2, T1 = 3,
                T2 = 0, T3 = 0, T4 = 0, L = 0, Inc = 0)
# run the particle filter with 20 particles
res2 <- c()
for(i in 1:10)
  res2[i] <- my_particleFilter(SEIT4L_stoch, theta, init.state,
                   data = FluTdC1971, n.particles = 35, brow = F)
plot(res, type = "l")
lines(res2, col = "red")


## Calibrate # of particles

# vector of number of particles to test
test.n.particles <- seq(50, 1000, 50)
# number of replicates 
n.replicates <- 100
# vector and data frame of results
sample.log.like <- vector("numeric", length = n.replicates)
res <- data.frame()

for(n.particles in test.n.particles){
  # start measuring time
  start.time  <- Sys.time()
  for(i in 1:n.replicates){
    cat("Testing ", i, " particles\n")
    # one Monte-Carlo estimate of the log-likelihood
    sample.log.like[i] <- my_particleFilter(SEIT4L_stoch, theta,
                                            init.state, FluTdC1971,
                                            n.particles)
  }
  # end measuring time
  end.time  <- Sys.time()
  
  # keep only replicate with finite log-likelihood to be able to compute the mean and sd
  # this give us the proportion of replicates with particle depletion.
  sample.finite.log.like <- sample.log.like[is.finite(sample.log.like)]
  
  ans <- c(mean = mean(sample.finite.log.like), 
           sd = sd(sample.finite.log.like), 
           prop.depletion = 1-length(sample.finite.log.like)/length(sample.log.like), 
           time = end.time - start.time)
  
  res <- rbind(res, t(ans))
}


## Setting pMCMC

# wrapper for posterior
my_posteriorSto <- function(theta){
  
  my_fitmodel <- SEIT4L_stoch
  my_init.state <- c("S" = 279,"E" = 0,"I" = 2,"T1" = 3,"T2" = 0, "T3" = 0, "T4" = 0,"L" = 0,"Inc" = 0) 
  
  my_n.particles <- 400 
  # you can reduce the number of particles if your pMCMC is too slow
  
  return(logPosterior(fitmodel = my_fitmodel,
                      theta = theta,
                      init.state = my_init.state,
                      data = FluTdC1971,
                      margLogLike = my_particleFilter,
                      n.particles = my_n.particles))
  
}

# load results of deterministic fit
data(mcmc_TdC_deter_longRun)

# Let's use the first trace only, no need to burn or thin
trace <- mcmc_SEITL_infoPrior_theta1$trace

# we will start the pMCMC at the mean posterior estimate
# of the deterministic fit
init.theta <- colMeans(trace[SEIT4L_stoch$theta.names])

# and we take the empirical covariance matrix for the 
# Gaussian kernel proposal
covmat <- mcmc_SEITL_infoPrior_theta1$covmat.empirical

# lower and upper limits of each parameter
lower <- c(R0 = 0, D_lat = 0 , D_inf = 0, alpha = 0, D_imm = 0, rho = 0)
upper <- c(R0 = Inf, D_lat = Inf , D_inf = Inf, alpha = 1, D_imm = Inf, rho = 1)

# number of iterations for the MCMC
n.iterations <- 50 # just a few since it takes quite a lot of time

# Here we don't adapt so that we can check the acceptance rate of the empirical covariance matrix
adapt.size.start <- 100
adapt.size.cooling <- 0.99
adapt.shape.start <- 100


# load traces
data(pmcmc_SEIT4L_infoPrior_n50)

# combine into a `mcmc.list` object
library("coda")
trace <- mcmc.list(lapply(pmcmc_SEIT4L_infoPrior_n50, function(chain) {
  mcmc(chain$trace)
}))

# acceptance rate is below the optimal 23%
1 - rejectionRate(trace)
# accordingly, the combined ESS is a bit low
effectiveSize(trace)
library("lattice")
xyplot(trace)
# this can take some time as we have 5 chains
plotESSBurn(trace)
# Actually, it looks like no burn-in is needed What about autocorrelation
acfplot(x = trace, lag.max = 50)

# There is substantial autocorrelation but we can't thin too much since the
# chains are quite short.  So let's keep 1 iteration every 20
trace.thin.n50 <- burnAndThin(trace, thin = 20)
# Finally we can plot the posterior density
densityplot(x = trace.thin.n50)

## 400 PArticles !!! 
data(pmcmc_SEIT4L_infoPrior_n400)
trace <- mcmc.list(lapply(pmcmc_SEIT4L_infoPrior_n400, function(chain) {
  mcmc(chain$trace)
}))
trace.thin.n400 <- burnAndThin(trace, thin = 20)
densityplot(x = trace.thin.n400)
plotPosteriorDensity(list(n50 = trace.thin.n50, n400 = trace.thin.n400))


## Stochastic vs deterministic fit 

# create mcmc object
trace1 <- mcmc(mcmc_SEIT4L_infoPrior_theta1$trace)
trace2 <- mcmc(mcmc_SEIT4L_infoPrior_theta2$trace)
# combine in a mcmc.list
trace <- mcmc.list(trace1, trace2)
# burn and thin as the chain with uniform prior (see above sections)
trace.deter <- burnAndThin(trace, burn = 5000, thin = 40)
# compare posterior density
plotPosteriorDensity(list(deter = trace.deter, sto = trace.thin.n400))

# combine all traces in a data frame
library("plyr")
trace.combined <- ldply(trace.thin.n400)
# take the mean of theta
theta.bar <- colMeans(trace.combined[SEIT4L_stoch$theta.names])
print(theta.bar)
# compute its log-likelihood
init.state <- c(S = 279, E = 0, I = 2, T1 = 3, T2 = 0, T3 = 0, T4 = 0, L = 0, 
                Inc = 0)
log.like.theta.bar <- my_particleFilter(SEIT4L_stoch, theta.bar, init.state, 
                                        data = FluTdC1971, n.particles = 40)
print(log.like.theta.bar)
log.like.theta.bar.deter <- dTrajObs(SEIT4L_deter, theta.bar, init.state,
                                     data = FluTdC1971, log = TRUE)
print(log.like.theta.bar.deter)
# and its deviance
D.theta.bar <- -2 * log.like.theta.bar
print(D.theta.bar)
# the effective number of parameters
p.D <- var(-2 * trace.combined$log.likelihood)/2
print(p.D)
# and finally the DIC
DIC <- D.theta.bar + 2 * p.D
print(DIC)

# take the mean posterior estimates of the deterministic model
x <- summary(trace.deter)
theta.bar.deter <- x$statistics[SEIT4L_deter$theta.names, "Mean"]

plotFit(SEIT4L_stoch, theta.bar, init.state, data = FluTdC1971, n.replicates = 1000)
plotFit(SEIT4L_deter, theta.bar.deter, init.state, data = FluTdC1971, n.replicates = 1000)



library("pomp")

density <- dwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
draw <- rwish(3, matrix(c(1,.3,.3,1),2,2))

plot(dwish( W = matrix(c(2,-.3,-.3,4),2,2), v = 5,
                      S = matrix(c(1,.3,.3,1),2,2)))
