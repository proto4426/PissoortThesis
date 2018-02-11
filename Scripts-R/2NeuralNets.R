load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1.Rdata")
#load("C:\\Users\\Piss\\Documents\\LINUX\\Documents\\Extreme\\R resources\\IRM\\data1.RData")

library(PissoortThesis)
library(GEVcdn)
library(gridExtra)
library(grid)
library(tidyverse)


## 1) Define the hierarchy of models of increasing complexity

models <- list()

## We follow the  approach of Martins and Stedinger (2000) for the prior parameters of
# the truncated beta density (beta.p = 9 and beta.q=6)
# Stationary model
models[[1]] <- list(Th = gevcdn.identity, beta.p = 9,  beta.q = 6,
                    fixed = c("location", "scale", "shape"))
# Linear models
models[[2]] <- list(Th = gevcdn.identity, beta.p = 9,  beta.q = 6,
                    fixed = c("shape","scale"))
models[[3]] <- list(Th = gevcdn.identity, beta.p = 9,  beta.q = 6,
                    fixed = c("shape"))
models[[4]] <- list(Th = gevcdn.identity, beta.p = 9,  beta.q = 6 )

# Nonlinear, 1 or 2 hidden nodes
models[[5]] <- list(n.hidden = 1, beta.p = 9,  beta.q = 6,
                    Th = gevcdn.logistic, fixed = c("shape", "scale"))
models[[6]] <- list(n.hidden = 2, beta.p = 9,  beta.q = 6,
                    Th = gevcdn.logistic, fixed = c("shape", "scale"))
models[[7]] <- list(n.hidden = 1, beta.p = 9,  beta.q = 6,
                    Th = gevcdn.logistic, fixed = c("shape"))
models[[8]] <- list(n.hidden = 2, beta.p = 9,  beta.q = 6,
                    Th = gevcdn.logistic, fixed = c("shape"))
models[[9]] <- list(n.hidden = 1, beta.p = 9,  beta.q = 6,
                    Th = gevcdn.logistic)
models[[10]] <- list(n.hidden = 2, beta.p = 9,  beta.q = 6,
                     Th = gevcdn.logistic)


## Try with the hyperbolic tangent activation function for nonlinear
hyp.tan <- function(x) ( exp(x) - exp(x) ) / ( exp(x) + exp(-x) )

models[[11]] <- list(n.hidden = 1, beta.p = 9,  beta.q = 6,
                     Th = hyp.tan, fixed = c("shape", "scale"))
models[[12]] <- list(n.hidden = 2, beta.p = 9,  beta.q = 6,
                     Th = hyp.tan, fixed = c("shape", "scale"))
models[[13]] <- list(n.hidden = 1, beta.p = 9,  beta.q = 6,
                     Th = hyp.tan, fixed = c("shape"))
models[[14]] <- list(n.hidden = 2, beta.p = 9,  beta.q = 6,
                     Th = hyp.tan, fixed = c("shape"))
models[[15]] <- list(n.hidden = 1, beta.p = 9,  beta.q = 6,
                     Th = hyp.tan)
models[[16]] <- list(n.hidden = 2, beta.p = 9,  beta.q = 6, Th = hyp.tan)




# Put the data to use in matrix x and y for ease of use inside GEVcdn framework
x <- as.matrix(seq(1, length(max_years$data)))
y <- as.matrix(max_years$data)


## 2) Fit the models and retrieve the weights
set.seed(123)
weights.models <- list()
for(i in seq_along(models)){
  weights.models[[i]] <- PissoortThesis::gevcdn.fit2(x = x, y = y,
                                                    n.trials = 1,
                                              n.hidden = models[[i]]$n.hidden,
                                              Th = models[[i]]$Th,
                                              fixed = models[[i]]$fixed )
}
# Printed outputs correspond to the value of the optimized ?gevcdn.cost
# (see also ?optim()),  for each model in this loop.


## Select best model

models.AICc <- round(sapply(weights.models, attr, which = "AICc"), 3)
# Comparing the AICc, we confirm that shape parameter must be held fixed.
# But last model seems also good...
models.BIC <- round(sapply(weights.models, attr, which = "BIC"), 3)
# Clear evidence for the 5th model (simple linear trend in location parameter)
# BIC penalizes more the more complex models --> pasimony
weights.best.aicc <- weights.models[[which.min(models.AICc)]]
weights.best.bic <- weights.models[[which.min(models.BIC)]]

parms.best <- gevcdn.evaluate(x, weights.best.bic)
# Find the parameter b_1 :
( parms.best[nrow(parms.best), "location"] - parms.best[1, "location"] ) / nrow(parms.best)


#### Compare nested models : Deviance statistics

## stationary vs linear trend
nll1 <- attr(weights.models[[1]], "NLL")
nll2 <- attr(weights.models[[2]], "NLL")
pchisq( 2 *( (-nll2) - (-nll1) ), lower.tail = F,
        df = attr(weights.models[[2]], "k") - attr(weights.models[[1]], "k"))
# ~exactly same result as previously done : linear trend clearly significant

# Compare the model chosen by AIC and by BIC.
nll9 <- attr(weights.models[[9]], "NLL")
pchisq( 2 *( (-nll9) - (-nll1) ), lower.tail = F,
        df = attr(weights.models[[9]], "k") - attr(weights.models[[2]], "k"))

## we can easily guess the output produced for the other tests



##### With weight regularization (way to prevent overfitting) :
weights.models_regul <- list()
for(i in seq_along(models)){
  weights.models_regul[[i]] <-
    PissoortThesis::gevcdn.fit2(x = x,
                                y = y,
                                n.trials = 1,
                                n.hidden = models[[i]]$n.hidden,
                                Th = models[[i]]$Th,
                                fixed = models[[i]]$fixed,
                                sd.norm = 5)
}
# By putting a prior on the weights, the number of effective parameters will change
# and hence, we cannot use AICc or BIC as above. The idea would be to select the
# best model and then use cross-validation to select optimal sd.norm. However,
# as our best model (BIC) has no hidden layers, it is not necessary to regularize.
# We will rather use bagging in the next step to introduce flexibility.


### Compute the quantiles an plot the results

'gev_cdnFitPlot' <- function(parms.best){
  
  x <- as.matrix(1901:2016) ; y <- as.matrix(max_years$data)
 
   q.best <- sapply(c(0.025, 0.05, 0.1,  0.5, 0.9, 0.95, 0.975), qgev,
                   location = parms.best[,"location"],
                   scale = parms.best[,"scale"],
                   shape = parms.best[,"shape"])
  
  
  # matplot(x, cbind(y, q.best), type = c("b", rep("l", 6)),
  #         lty = c(1, rep(c(1, 2, 1), 2)),
  #         lwd = c(1, rep(c(3, 2, 3), 2)),
  #         col = c("black", rep("blue", 3), rep("blue", 3)),
  #         pch = 19, xlab = "x", ylab = "y", main = "gevcdn.fit")
  
  # Or in ggplot
  df <- data.frame(year = x , obs = y, q.025 = q.best[,1], q.05 = q.best[,2],
                   q.10 = q.best[,3],
                   q.50 = q.best[,4], q.90 = q.best[,5], q.975 = q.best[,7])
  col.quantiles <- c("2.5% and 97.5%" = "blue", "10% and 90%" = "green", "50%" = "red")
  gg.cdn <- ggplot(df, aes(x = year, y = obs)) +
    geom_line() + geom_point() +
    geom_line(aes(y = q.025, col = "2.5% and 97.5%")) +
    geom_line(aes(y = q.50, col = "50%")) +
    geom_line(aes(y = q.975, col = "2.5% and 97.5%")) +
    geom_line(aes(y = q.10, col = "10% and 90%")) +
    geom_line(aes(y = q.90, col = "10% and 90%")) +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    labs(title = "GEV-CDN quantiles with identity link for the location µ(t)",
         y =  expression( Max~(T~degree*C))) +
    theme_piss()
  return(list(gg.cdn, df))
}

gg.cdn <- gev_cdnFitPlot(parms.best = parms.best)[[1]]


## For the 2 hidden layers model (model) comparison
parms.best7 <- gevcdn.evaluate(x, weights.models[[7]])

ggcdn7 <- gev_cdnFitPlot(parms.best = parms.best7)


#####################################
#### BAGGING #######################
####################################


n.boot <- 10 # Number of boostrapped iterations

time <- proc.time()
weights.on <- gevcdn.bag(x = x, y = y, iter.max = 100,
                         iter.step = 10, n.bootstrap = n.boot,
                        n.hidden = 2)
(proc.time()-time)[3]   # 8.4 sec on our machine for n.boot <- 10

parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)


## Create a function to calculate average for the parameters accross the bootstrap resamples (in list param)
'mean_of_list' <- function(param = parms.on, brow = F) {
  if(brow) browser()
  parms <- matrix(0, nrow = nrow(param[[1]]), ncol = 3)
  for (i in 1:length(param)){
    parms <- parms + as.matrix(param[[i]])
  }
  parms <- parms / length(param)
  parms
  #parms <- apply(parms.on, 2, FUN = mean)  # Try with this : vectorized version
  #colMeans(do.call(rbind, parms.on))
}


## Do it through parallel computing since it is very slow if we want a correct nuumber of resamp
library(foreach)
library(doParallel)


'bag_ggplot' <- function(M = 50, n.boot = 10){
  
  x <- as.matrix(1901:2016)
  y <- as.matrix(max_years$data)
  
  #setup parallel backend to use many processors
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
  registerDoParallel(cl)
  
  #M <- 50 # M is the number of times we do bagging, and n.boot is the number
  #of resamples inside gevcdn.bag()
  'mean_of_list' <- function(param = parms.on, brow = F) {
    if(brow) browser()
    parms <- matrix(0, nrow = nrow(param[[1]]), ncol = 3)
    for (i in 1:length(param)){
      parms <- parms + as.matrix(param[[i]])
    }
    parms <- parms / length(param)
    parms
  }
  t <- proc.time()
  bag_par <- foreach(i = 1:M,
                     .packages = c("PissoortThesis", "GEVcdn"),
                     .verbose = T) %dopar% {
                       set.seed(i+1234)
                       weights.on <- gevcdn.bag(x = x, y = y,
                                                iter.max = 100,
                                                fixed = c("shape"),
                                                iter.step = 10,
                                                n.bootstrap = n.boot,
                                                n.hidden = 2,
                                                sd.norm = 5)
                       parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)
                       
                       mean_of_list(parms.on)
                     }
  (proc.time()-t)[3]  # 17sec for 5*10=50 resamp. // 158sec for 50*10=500 resamp.
  stopCluster(cl)
  
  bag_par_on <- mean_of_list(bag_par)
  
  
  ## Estimate the location parameter b_1
  ( bag_par_on[nrow(bag_par_on), "location"] - bag_par_on[1, "location"] ) / nrow(bag_par_on)
  ## And parameter alpha_0 and alpha_1
  alpha_0 <- log(bag_par_on[1,2])   ;   alpha_1 <- c()
  for (i in 2:(nrow(bag_par_on)) ) {
    alpha_1[i] <-  ( log( bag_par_on[i, "scale"] ) - alpha_0 ) / i
  }
  
  
  q <- t(sapply(max_years$data, quantile, probs = c(.1, .5, .9)))
  
  q.025.on <- q.05.on <- q.10.on <- q.50.on <- q.90.on <- q.95.on <- q.975.on <- c()
  for(i in seq_along(bag_par)){
    q.025.on <- cbind(q.025.on, VGAM::qgev(p = 0.025,
                                           location = bag_par[[i]][,"location"],
                                           scale = bag_par[[i]][,"scale"],
                                           shape = bag_par[[i]][,"shape"]))
    q.05.on <- cbind(q.05.on, VGAM::qgev(p = 0.05,
                                         location = bag_par[[i]][,"location"],
                                         scale = bag_par[[i]][,"scale"],
                                         shape = bag_par[[i]][,"shape"]))
    q.10.on <- cbind(q.10.on, VGAM::qgev(p = 0.1,
                                         location = bag_par[[i]][,"location"],
                                         scale = bag_par[[i]][,"scale"],
                                         shape = bag_par[[i]][,"shape"]))
    q.50.on <- cbind(q.50.on, VGAM::qgev(p = 0.5,
                                         location = bag_par[[i]][,"location"],
                                         scale = bag_par[[i]][,"scale"],
                                         shape = bag_par[[i]][,"shape"]))
    q.90.on <- cbind(q.90.on, VGAM::qgev(p = 0.9,
                                         location = bag_par[[i]][,"location"],
                                         scale = bag_par[[i]][,"scale"],
                                         shape = bag_par[[i]][,"shape"]))
    q.95.on <- cbind(q.975.on, VGAM::qgev(p = 0.95,
                                          location = bag_par[[i]][,"location"],
                                          scale = bag_par[[i]][,"scale"],
                                          shape = bag_par[[i]][,"shape"]))
    q.975.on <- cbind(q.975.on, VGAM::qgev(p = 0.975,
                                           location = bag_par[[i]][,"location"],
                                           scale = bag_par[[i]][,"scale"],
                                           shape = bag_par[[i]][,"shape"]))
  }
  
  ## Plot data and quantiles
  # matplot(cbind(y, q, rowMeans(q.05.on), rowMeans(q.10.on), rowMeans(q.50.on),
  #               rowMeans(q.90.on), rowMeans(q.95.on)), type = c("b", rep("l", 10)),
  #         lty = c(1, rep(c(1, 2, 1), 4)),
  #         lwd = c(1, rep(c(3, 2, 3), 4)),
  #         col = c("red", rep("orange", 3), "blue", "green", "black", "green","blue"),
  #         pch = 19, xlab = "x", ylab = "y",
  #         main = "gevcdn.bag (early stopping on)")
  
  # With ggplot
  df.bag <- data.frame(year = x , obs = y, q.025 = rowMeans(q.025.on),
                       q.10 = rowMeans(q.10.on), q.50 = rowMeans(q.50.on),
                       q.90 = rowMeans(q.90.on), q.975 = rowMeans(q.975.on))
  col.quantiles <- c("2.5% and 97.5%" = "blue", "10% and 90%" = "green", "50%" = "red")
  gg.bag <- ggplot(df.bag, aes(x = year, y = obs)) +
    geom_line() + geom_point() +
    geom_line(aes(y = q.025, col = "2.5% and 97.5%")) +
    geom_line(aes(y = q.50, col = "50%")) +
    geom_line(aes(y = q.975, col = "2.5% and 97.5%")) +
    geom_line(aes(y = q.10, col = "10% and 90%")) +
    geom_line(aes(y = q.90, col = "10% and 90%")) +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    labs(title = "GEV-CDN bagging (early stopping, 2 hidden nodes) quantiles", y = "") +
    theme_piss()
  gg.bag
}

gg.bag <- bag_ggplot()
gg.bag5 <- bag_ggplot(1, 2)


## Gather Nostationary 2 hidden layers model with and without bagging in one plot. 
## See the over-complexity when no bagging step is made.

gg.bag2 <- gg.bag + geom_line(data = ggcdn7[[2]], aes(x = year, y = q.025), col = "blue", linetype = "dashed") + 
  geom_line(data = ggcdn7[[2]], aes(x = year, y = q.50), col = "red", linetype = "dashed") +
  geom_line(data = ggcdn7[[2]], aes(x = year, y = q.975), col = "blue", linetype = "dashed") +
  geom_line(data = ggcdn7[[2]], aes(x = year, y = q.10), col = "green", linetype = "dashed") +
  geom_line(data = ggcdn7[[2]], aes(x = year, y = q.90), col = "green", linetype = "dashed") +
  # scale_x_continuous(breaks = c(seq(1900,2000, by = 20)), 
  #                    label = c(seq(1900,2000, by = 20))) + 
  labs(title = "GEV-CDN (bagging) quantiles for model with 2 hidden nodes ", y = "") +
       #subtitle = "Dashed lines represent the quantiles of the model with NO bagging : tendency to overfit is clear") +
  theme(plot.subtitle = element_text(18, hjust = 0.5) )



## Together the GEV-CDN without and with bagging
gridExtra::grid.arrange(gg.cdn, gg.bag, nrow = 1)
PissoortThesis::grid_arrange_legend(gg.cdn, gg.bag2)

## Comparison of the quantiles in one Figure :

df_lin_bag <- cbind.data.frame(df = df, df.bag = df.bag, Year = 1901:2016)

ggplot(df_lin_bag, aes(x = Year)) +
  geom_line(aes(y = df.obs), size = 0.2) + geom_point(aes(y = df.obs), size = 0.5) +
  geom_line(aes(y = df.q.025, col = "2.5% and 97.5%")) +
  geom_line(aes(y = df.q.50, col = "50%")) +
  geom_line(aes(y = df.q.975, col = "2.5% and 97.5%")) +
  geom_line(aes(y = df.q.10, col = "10% and 90%")) +
  geom_line(aes(y = df.q.90, col = "10% and 90%")) +
  geom_line(aes(y = df.bag.q.025, col = "2.5% and 97.5%"), linetype = "dashed") +
  geom_line(aes(y = df.bag.q.50, col = "50%"), linetype = "dashed") +
  geom_line(aes(y = df.bag.q.975, col = "2.5% and 97.5%"), linetype = "dashed") +
  geom_line(aes(y = df.bag.q.10, col = "10% and 90%"), linetype = "dashed") +
  geom_line(aes(y = df.bag.q.90, col = "10% and 90%"), linetype = "dashed") +
  scale_colour_manual(name = "Quantiles", values = col.quantiles) +
  labs(title = "Annual maxima and Quantiles of a linear and a nonlinear GEV-CDN model with Bagging",
       y = expression( Max~(T~degree*C)),
       subtitle = "Lines represent quantiles of the linear model on the location and dashed represent quantiles of the bagged nonlinear model on location and scale") +
  theme_piss(legend.position = c(.947, .105)) +
  theme(plot.subtitle = text(18, hjust = 0.5) )


## Comparisons of the fitted quantiles with the empirical quantiles
quantile(max_years$df$Max, probs = c(0.05, 0.1, 0.5, 0.9, 0.95))
q.best   ;    df.bag [,3:7]




###### Bootstrap confidence intervals ########
#############################################



## Fit 30 bootstrapped models
time <- proc.time()
set.seed((1234))
ci.lin <- gevcdn.bootstrap(n.bootstrap = 500,
                           x = x, y = y,
                         iter.max = 100,
                         n.hidden = 0,
                         fixed = c("shape", "scale"),
                         Th = gevcdn.identity,
                         n.trials = 1,
                         boot.method = "residual",
                         probs = c(0.1, 0.5, 0.9))
(proc.time() - time)[3]  # 11sec for B=50 , 61sec for B=250





# M <- 100 # M is the number of times we do bagging, thus the total number of resamples is
# n.boot*M.
# Function created just to facilitate our work to compute directly several methods.
# LAter include inside the package.
"super_boot.ci_parallel" <- function(M = 100, boot_method = "residual"){

  x <- as.matrix(seq(1, length(max_years$data)))
  y <- as.matrix(max_years$data)

  ### In parralel !!! First model : nonstationary in location ( 0 hidden)
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
  doParallel::registerDoParallel(cl)

  n.boot <- 10  #, and n.boot is the number of resamples inside gevcdn.bag()
  t <- proc.time()
  boot_par <- foreach::foreach(i = 1:M,
                              .packages = c("PissoortThesis", "GEVcdn"),
                              .verbose = T) %dopar% {
                        set.seed(i+12)
                        ci.lin <- gevcdn.bootstrap(n.bootstrap = n.boot,
                                                   x = x, y = y,
                                                   iter.max = 100,
                                                   n.hidden = 0,
                                                   fixed = c("shape", "scale"),
                                                   Th = gevcdn.identity,
                                                   n.trials = 1,
                                                   boot.method = boot_method,
                                                   probs = c(0.1, 0.5, 0.9))
                        list(loc = ci.lin$location.bootstrap,
                             sc = ci.lin$scale.bootstrap,
                             sh = ci.lin$shape.bootstrap
                        )
                      }
  (proc.time()-t)[3] ## only 74 sec for B=1000 !!!
  parallel::stopCluster(cl)

  ## Aggregate the results of the resamples and the parralel computation

  boot.loc <- matrix(nrow = nrow(boot_par[[1]]$loc))
  for (i in 1:M)   boot.loc <- cbind(boot.loc, boot_par[[i]]$loc)
  boot.loc_f <- boot.loc[,-1]

  boot.sc <- matrix(nrow = nrow(boot_par[[1]]$loc))
  for (i in 1:M)   boot.sc <- cbind(boot.sc, boot_par[[i]]$sc)
  boot.sc_f <- boot.sc[,-1]

  boot.sh <- matrix(nrow = nrow(boot_par[[1]]$loc))
  for (i in 1:M)   boot.sh <- cbind(boot.sh, boot_par[[i]]$sh)
  boot.sh_f <- boot.sh[,-1]

  ## Compute the quantiles
  b.loc <- t(apply(boot.loc_f, 1, quantile, p = c(0.025, 0.975)))
  b.sc <- t(apply(boot.sc_f, 1, quantile, p = c(0.025, 0.975)))
  b.sh <- t(apply(boot.sh_f, 1, quantile, p = c(0.025, 0.975)))


  #### Complex model  : 2 hidden layers for location and scale
  ###################

  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
  doParallel::registerDoParallel(cl)

  n.boot <- 10
  # M <- 100 # M is the number of times we do bagging, and n.boot is the number
  #of resamples inside gevcdn.bag()
  t <- proc.time()
  boot_par_2hidden <- foreach::foreach(i = 1:M,
                                      .packages = c("PissoortThesis", "GEVcdn"),
                                      .verbose = T) %dopar% {
                        set.seed(i+12)
                                        ci.lin <- gevcdn.bootstrap(n.bootstrap = n.boot,
                                                                   x = x, y = y,
                                                                   iter.max = 100,
                                                                   n.hidden = 2,
                                                                   fixed = c("shape"),
                                                                   Th = gevcdn.logistic,
                                                                   n.trials = 1,
                                                                   boot.method = boot_method,
                                                                   probs = c(0.1, 0.5, 0.9))
                                        list(loc = ci.lin$location.bootstrap,
                                             sc = ci.lin$scale.bootstrap,
                                             sh = ci.lin$shape.bootstrap
                                        )

                              }
  (proc.time()-t)[3] ## only 74 sec for B=1000 !!!
  parallel::stopCluster(cl)

  ## Aggregate the results of the resamples and the parralel computation

  boot.loc.h <- matrix(nrow = nrow(boot_par_2hidden[[1]]$loc))
  for (i in 1:M)   boot.loc.h <- cbind(boot.loc.h, boot_par_2hidden[[i]]$loc)
  boot.loc.h_f <- boot.loc.h[,-1]

  boot.sc.h <- matrix(nrow = nrow(boot_par_2hidden[[1]]$loc))
  for (i in 1:M)   boot.sc.h <- cbind(boot.sc.h, boot_par_2hidden[[i]]$sc)
  boot.sc.h_f <- boot.sc.h[,-1]

  boot.sh.h <- matrix(nrow = nrow(boot_par_2hidden[[1]]$loc))
  for (i in 1:M)   boot.sh.h <- cbind(boot.sh.h, boot_par_2hidden[[i]]$sh)
  boot.sh.h_f <- boot.sh.h[,-1]

  ## Compute the quantiles
  b.loc.h <- t(apply(boot.loc.h_f, 1, quantile, p = c(0.025, 0.975)))
  b.sc.h <- t(apply(boot.sc.h_f, 1, quantile, p = c(0.025, 0.975)))
  b.sh.h <- t(apply(boot.sh.h_f, 1, quantile, p = c(0.025, 0.975)))

  # Data frames for the two consider models, simple and complex(.h)
  df.boot <- data.frame(year = 1901:2016 , obs = y,
                        loc = b.loc,
                        scale = b.sc,
                        shape = b.sh)
  df.boot.h <- data.frame(year = 1901:2016 , obs = y,
                          loc = b.loc.h,
                          scale = b.sc.h,
                          shape = b.sh.h)


  ### FIRST model :
  ## Plot data and percentile confidence intervals for GEV parameters

  # for the location
  gg.boot.loc <- ggplot(df.boot, aes(x = year)) +
    geom_line(aes(y = loc.2.5.), col = "blue") +
    geom_line(aes(y = loc.97.5.), col = "blue") +
    labs(title = "GEV-CDN identity link on location only (0 hidden)", y = "") +
    theme_piss() +
    coord_cartesian(ylim = c(min(df.boot$loc.2.5., df.boot.h$loc.2.5.), # Put results on same scale
                             max(df.boot.h$loc.97.5., df.boot$loc.97.5.)))
  # for the scale
  gg.boot.sc <- ggplot(df.boot, aes(x = year)) +
    geom_line(aes(y = scale.2.5.), col = "green") +
    geom_line(aes(y = scale.97.5.), col = "green") +
    labs(title = "GEV-CDN identity link on location only (0 hidden)", y = "") +
    theme_piss() +
    coord_cartesian(ylim = c(min(df.boot$scale.2.5., df.boot.h$scale.2.5.),  # put results on same scale
                             max(df.boot.h$scale.97.5., df.boot$scale.97.5.)))
  # for the shape
  gg.boot.sh <- ggplot(df.boot, aes(x = year)) +
    geom_line(aes(y = shape.2.5.), col = "red") +
    geom_line(aes(y = shape.97.5.), col = "red") +
    labs(title = "GEV-CDN identity link on location only (0 hidden)", y = "") +
    theme_piss() +
    coord_cartesian(ylim = c(min(df.boot$shape.2.5., df.boot.h$shape.2.5.), # Put results on same scale
                             max(df.boot.h$shape.97.5., df.boot$shape.97.5.)))

  # All in one for the identity link model on the location
  g_lin_boot <- grid.arrange(gg.boot.loc, gg.boot.sc, gg.boot.sh, nrow = 2,
                             top = textGrob(expression("GEV-CDN linear on the location only (0 hidden)"),
                                            gp = gpar(fontsize = 17, font = 4, col ="black")))
  g_lin_boot



  ### COMPLEX (second) model
  ## Plot data and percentile confidence intervals for GEV parameters

  # For the location
  gg.boot.loc.h <- ggplot(df.boot.h, aes(x = year)) +
    geom_line(aes(y = loc.2.5.), col = "blue") +
    geom_line(aes(y = loc.97.5.), col = "blue") +
    labs(title = "GEV-CDN nonlinear sigmoid with 2 hidden layers ", y = "") +
    theme_piss() +
    coord_cartesian(ylim = c(min(df.boot$loc.2.5., df.boot.h$loc.2.5.),  # Put results on same scale
                             max(df.boot.h$loc.97.5., df.boot$loc.97.5.)))
  # For the scale
  gg.boot.sc.h <- ggplot(df.boot.h, aes(x = year)) +
    geom_line(aes(y = scale.2.5.), col = "green") +
    geom_line(aes(y = scale.97.5.), col = "green") +
    labs(title = "GEV-CDN nonlinear sigmoid with 2 hidden layers ", y = "") +
    theme_piss() +
    coord_cartesian(ylim = c(min(df.boot$scale.2.5., df.boot.h$scale.2.5.),  # Put results on same scale
                             max(df.boot.h$scale.97.5., df.boot$scale.97.5.)))
  # For the shape
  gg.boot.sh.h <- ggplot(df.boot.h, aes(x = year)) +
    geom_line(aes(y = shape.2.5.), col = "red") +
    geom_line(aes(y = shape.97.5.), col = "red") +
    labs(title = "GEV-CDN nonlinear sigmoid with 2 hidden layers ", y = "") +
    theme_piss() +
    coord_cartesian(ylim = c(min(df.boot$shape.2.5., df.boot$shape.2.5.),  # Put results on same scale
                             max(df.boot.h$shape.97.5., df.boot$shape.97.5.)))


  ## All in one for the complex model
  g_hidd <- grid.arrange(gg.boot.loc.h, gg.boot.sc.h, gg.boot.sh.h, nrow = 2,
                         top = textGrob(expression("GEV-CDN nonlinear on location and scale (2 hidden)"),
                                        gp = gpar(fontsize = 17, font = 4, col ="black")))
  g_hidd


  ## ALL together
  grid.arrange(g_lin_boot, g_hidd)

  ## Return a list to be able to plot the confidence intervals by parameters
  return(list(g.loc = gg.boot.loc, g.loc.h = gg.boot.loc.h,
              g.sc =  gg.boot.sc, g.sc.h = gg.boot.sc.h,
              g.sh = gg.boot.sh, g.sh.h = gg.boot.sh.h,

              df.lin = df.boot, df.hidden = df.boot.h))
}



## All results for the RESIDUAL bootstrap
res.boot.residual <- super_boot.ci_parallel(M = 10, boot_method = "residual")

grid_loc_res <-  grid.arrange(res.boot.residual$g.loc, res.boot.residual$g.loc.h, nrow = 1,
             top = textGrob(expression(" 95% residual bootstrap confidence interval for the Location parameter"),
                            gp = gpar(fontsize = 19, font = 2, col ="black")))
grid_sc_res <- grid.arrange(res.boot.residual$g.sc, res.boot.residual$g.sc.h, nrow = 1,
             top = textGrob(expression(" 95% residual bootstrap confidence interval for the Scale parameter"),
                            gp = gpar(fontsize = 19, font = 2, col ="black")))
grid_sh_res <- grid.arrange(res.boot.residual$g.sh, res.boot.residual$g.sh.h, nrow = 1,
             top = textGrob(expression("  95% residual bootstrap confidence interval for the Shape parameter"),
                            gp = gpar(fontsize = 19, font = 2, col ="black")))
grid.arrange(grid_loc_res, grid_sc_res, grid_sh_res)


## And for the PARAMETRIC bootstrap
res.boot.par <- super_boot.ci_parallel(M = 10, boot_method = "parametric")

grid_loc_par <-  grid.arrange(res.boot.par$g.loc, res.boot.par$g.loc.h, nrow = 1,
                              top = textGrob(expression(" 95% parametric bootstrap confidence interval for the Location parameter"),
                                             gp = gpar(fontsize = 19, font = 2, col ="black")))
grid_sc_par <- grid.arrange(res.boot.par$g.sc, res.boot.par$g.sc.h, nrow = 1,
                            top = textGrob(expression(" 95% parametric bootstrap confidence interval for the Scale parameter"),
                                           gp = gpar(fontsize = 19, font = 2, col ="black")))
grid_sh_par <- grid.arrange(res.boot.par$g.sh, res.boot.par$g.sh.h, nrow = 1,
                            top = textGrob(expression(" 95% parametric bootstrap confidence interval for the Shape parameter"),
                                           gp = gpar(fontsize = 19, font = 2, col ="black")))
grid.arrange(grid_loc_par, grid_sc_par, grid_sh_par)



## Quantitative comparisons between the residual and parametric bootstraps

diff.boot.met.lin <- res.boot.residual$df.lin[,-(1:2)] - res.boot.par$df.lin[, -(1:2)]
diff.boot.met.hid <- res.boot.residual$df.hidden[,-(1:2)] - res.boot.par$df.hidden[, -(1:2)]


# Create a function to easily compare all the methods. df contains the residuals for all quantiles
# To include in the package
'compare_methods_quantiles_plot' <- function(df) {

  col.quantiles <- c("97.5%" = "red", "2.5%" = "blue")

  gloc.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
    geom_line(aes(y = loc.2.5., col  = "2.5%")) +
    geom_line(aes(y = loc.97.5., col  = "97.5%"))+
    labs(title = "Location parameter") +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    theme_piss()

  gsc.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
    geom_line(aes(y = scale.2.5., col  = "2.5%")) +
    geom_line(aes(y = scale.97.5., col  = "97.5%")) +
    labs(title = "Scale parameter") +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    theme_piss()

  gsh.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
    geom_line(aes(y = shape.2.5., col  = "2.5%")) +
    geom_line(aes(y = shape.97.5., col  = "97.5%")) +
    labs(title = "Shape parameter") +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    theme_piss()
  #grid.arrange(gloc.lin, gsc.lin, gsh.lin, nrow = 1)
  PissoortThesis::grid_arrange_legend(gloc.lin, gsc.lin, gsh.lin)
}

# For the linear model : residual vs parametric bootstrap, for each parameters
compare_methods_quantiles_plot(diff.boot.met.lin)
# For the complex model : residual vs parametric bootstrap, for each parameters
compare_methods_quantiles_plot(diff.boot.met.hid)



### ==================================================================================
## Find the values of the confidence intervals for comparisons in Bayes_own_gev.R

cores <- parallel::detectCores()  ;   cl <- parallel::makeCluster(cores[1]-1)
doParallel::registerDoParallel(cl)  ; n.boot <- 10  ;  M <- 50  ;  t <- proc.time()
boot_par0 <- foreach::foreach(i = 1:M,
                             .packages = c("PissoortThesis", "GEVcdn"),
                             .verbose = T) %dopar% {
                               set.seed(i+12)
                               ci.lin <- gevcdn.bootstrap(n.bootstrap = n.boot,
                                                          x = x, y = y,
                                                          iter.max = 100,
                                                          n.hidden = 0,
                                                          fixed = c("shape", "scale"),
                                                          Th = gevcdn.identity,
                                                          n.trials = 1,
                                                          boot.method = "residual",
                                                          probs = c(0.1, 0.5, 0.9))
                               list(loc = ci.lin$location.bootstrap,
                                    sc = ci.lin$scale.bootstrap,
                                    sh = ci.lin$shape.bootstrap
                               )
                             }
(proc.time()-t)[3] ## only 74 sec for B=1000 !!!
parallel::stopCluster(cl)

## Aggregate the results of the resamples and the parralel computation
boot.loc <- matrix(nrow = nrow(boot_par[[1]]$loc))
for (i in 1:M)   boot.loc <- cbind(boot.loc, boot_par0[[i]]$loc)
boot.loc_f <- boot.loc[,-1]

boot.sc <- matrix(nrow = nrow(boot_par[[1]]$loc))
for (i in 1:M)   boot.sc <- cbind(boot.sc, boot_par0[[i]]$sc)
boot.sc_f <- boot.sc[,-1]

boot.sh <- matrix(nrow = nrow(boot_par[[1]]$loc))
for (i in 1:M)   boot.sh <- cbind(boot.sh, boot_par0[[i]]$sh)
boot.sh_f <- boot.sh[,-1]

## Compute the quantiles
t(apply(boot.loc_f, 1, quantile, p = c(0.025, 0.975)))
# For the location, as it takes the aggregated parameter mu(t), we cannot compute one single interval
sigma95_boot <- t(apply(boot.sc_f, 1, quantile, p = c(0.025, 0.5, 0.975)))[1,]
xi95_boot <- t(apply(boot.sh_f, 1, quantile, p = c(0.025, 0.5, 0.975)))[1,]
sigma75_boot <- t(apply(boot.sc_f, 1, quantile, p = c(0.25, 0.5,  0.75)))[1,]
xi75_boot <- t(apply(boot.sh_f, 1, quantile, p = c(0.25, 0.5, 0.75)))[1,]

boot.ci.Sig_Xi <- data.frame(sigma95_boot, sigma75_boot,
                             xi75_boot, xi95_boot)



## Find the values of the confidence intervals for comparisons in Bayes_own_gev.R
## Parametric bootstrap
cores <- parallel::detectCores()  ;
cl <- parallel::makeCluster(cores[1]-1)
doParallel::registerDoParallel(cl)  ; n.boot <- 10  ;  M <- 50  ;  t <- proc.time()
boot_par <- foreach::foreach(i = 1:M,
                             .packages = c("PissoortThesis", "GEVcdn"),
                             .verbose = T) %dopar% {
                               set.seed(i+12)
                               ci.lin <- gevcdn.bootstrap(n.bootstrap = n.boot,
                                                          x = x, y = y,
                                                          iter.max = 100,
                                                          n.hidden = 0,
                                                          fixed = c("shape", "scale"),
                                                          Th = gevcdn.identity,
                                                          n.trials = 1,
                                                          boot.method = "parametric",
                                                          probs = c(0.1, 0.5, 0.9))
                               list(loc = ci.lin$location.bootstrap,
                                    sc = ci.lin$scale.bootstrap,
                                    sh = ci.lin$shape.bootstrap
                               )
                             }
(proc.time()-t)[3] ## only 74 sec for B=1000 !!!
parallel::stopCluster(cl)

## Aggregate the results of the resamples and the parralel computation
boot.loc <- matrix(nrow = nrow(boot_par[[1]]$loc))
for (i in 1:M)   boot.loc <- cbind(boot.loc, boot_par[[i]]$loc)
boot.loc_f <- boot.loc[,-1]

boot.sc <- matrix(nrow = nrow(boot_par[[1]]$loc))
for (i in 1:M)   boot.sc <- cbind(boot.sc, boot_par[[i]]$sc)
boot.sc_f <- boot.sc[,-1]

boot.sh <- matrix(nrow = nrow(boot_par[[1]]$loc))
for (i in 1:M)   boot.sh <- cbind(boot.sh, boot_par[[i]]$sh)
boot.sh_f <- boot.sh[,-1]

## Compute the quantiles
t(apply(boot.loc_f, 1, quantile, p = c(0.025, 0.975)))
# For the location, as it takes the aggregated parameter mu(t), we cannot compute one single interval
sigma95_bootp <- t(apply(boot.sc_f, 1, quantile, p = c(0.025, 0.5, 0.975)))[1,]
xi95_bootp <- t(apply(boot.sh_f, 1, quantile, p = c(0.025, 0.5, 0.975)))[1,]
sigma75_bootp <- t(apply(boot.sc_f, 1, quantile, p = c(0.25, 0.5,  0.75)))[1,]
xi75_bootp <- t(apply(boot.sh_f, 1, quantile, p = c(0.25, 0.5, 0.75)))[1,]

boot.ci.Sig_Xi.par <- data.frame(sigma95_bootp, sigma75_bootp,
                             xi75_bootp, xi95_bootp)

