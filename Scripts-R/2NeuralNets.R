load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1.Rdata")
#load("C:\\Users\\Piss\\Documents\\LINUX\\Documents\\Extreme\\R resources\\IRM\\data1.RData")

library(PissoortThesis)
library(GEVcdn)


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
set.seed(1234)
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

q.best <- sapply(c(0.05, 0.1,  0.5, 0.9, 0.95), qgev,
                 location = parms.best[,"location"],
                 scale = parms.best[,"scale"],
                 shape = parms.best[,"shape"])


matplot(x, cbind(y, q.best), type = c("b", rep("l", 6)),
        lty = c(1, rep(c(1, 2, 1), 2)),
        lwd = c(1, rep(c(3, 2, 3), 2)),
        col = c("black", rep("blue", 3), rep("blue", 3)),
        pch = 19, xlab = "x", ylab = "y", main = "gevcdn.fit")

# Or in ggplot
df <- data.frame(year = x , obs = y, q.05 = q.best[,1], q.10 = q.best[,2],
                 q.50 = q.best[,3], q.90 = q.best[,4], q.95 = q.best[,5])
gg.cdn <- ggplot(df, aes(x = year, y = obs)) +
  geom_line() + geom_point() +
  geom_line(aes(y = q.05), col = "blue") +
  geom_line(aes(y = q.50), col = "red", linetype = "dashed") +
  geom_line(aes(y = q.95), col = "blue") +
  geom_line(aes(y = q.10), col = "green") +
  geom_line(aes(y = q.90), col = "green") +
  labs(title = "GEV-CDN quantiles with identity link for the location µ(t)",
       y =  expression(Centered ~ Max~(T~degree*C))) +
  theme_piss()
gg.cdn




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

#setup parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
registerDoParallel(cl)

M <- 2 # M is the number of times we do bagging, and n.boot is the number
#of resamples inside gevcdn.bag()
t <- proc.time()
bag_par <- foreach(i = 1:M,
                   .packages = c("PissoortThesis", "GEVcdn"),
                   .verbose = T) %dopar% {
           weights.on <- gevcdn.bag(x = x, y = y,
                                    iter.max = 100,
                                    fixed = c("shape"),
                                    iter.step = 10,
                                    n.bootstrap = n.boot,
                                    n.hidden = 2,
                                    sd.norm = .2)
           parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)

           mean_of_list(parms.on)
    }
(proc.time()-t)[3]  # 17sec for 5*10=50 resamp. // 158sec for 50*10=500 resamp.
stopCluster(cl)

bag_par_on <- mean_of_list(bag_par)


q <- t(sapply(max_years$data, quantile, probs = c(.1, .5, .9)))

q.05.on <- q.10.on <- q.50.on <- q.90.on <- q.95.on <- c()
for(i in seq_along(bag_par)){
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
  q.95.on <- cbind(q.95.on, VGAM::qgev(p = 0.95,
                                 location = bag_par[[i]][,"location"],
                                 scale = bag_par[[i]][,"scale"],
                                 shape = bag_par[[i]][,"shape"]))
}

## Plot data and quantiles
matplot(cbind(y, q, rowMeans(q.05.on), rowMeans(q.10.on), rowMeans(q.50.on),
              rowMeans(q.90.on), rowMeans(q.95.on)), type = c("b", rep("l", 10)),
        lty = c(1, rep(c(1, 2, 1), 4)),
        lwd = c(1, rep(c(3, 2, 3), 4)),
        col = c("red", rep("orange", 3), "blue", "green", "black", "green","blue"),
        pch = 19, xlab = "x", ylab = "y",
        main = "gevcdn.bag (early stopping on)")

# With ggplot
df.bag <- data.frame(year = x , obs = y, q.05 = rowMeans(q.05.on),
                     q.10 = rowMeans(q.10.on), q.50 = rowMeans(q.50.on),
                     q.90 = rowMeans(q.90.on), q.95 = rowMeans(q.95.on))
gg.bag <- ggplot(df.bag, aes(x = year, y = obs)) +
  geom_line() + geom_point() +
  geom_line(aes(y = q.05), col = "blue") +
  geom_line(aes(y = q.50), col = "red", linetype = "dashed") +
  geom_line(aes(y = q.95), col = "blue") +
  geom_line(aes(y = q.10), col = "green") +
  geom_line(aes(y = q.90), col = "green") +
  labs(title = "GEV-CDN bagging (early stopping, 2 hidden nodes) quantiles",
       y = "") +
  theme_piss()
gg.bag


## Together the GEV-CDN without and with bagging
gridExtra::grid.arrange(gg.cdn, gg.bag, nrow = 1)




###### Bootstrap confidence intervals ########

## Fit 30 bootstrapped models
CI <- gevcdn.bootstrap(n.bootstrap = 30, x = x, y = y,
                       iter.max = 100, n.hidden = 2,
                       Th = gevcdn.logistic, n.trials = 1,
                       boot.method = "residual",
                       probs = c(0.1, 0.5, 0.9))
## Plot data and percentile confidence intervals for GEV parameters
par(mfrow = c(2, 2))
matplot(x, y, type = "b", pch = 19, col = "red", xlab = "x",
        ylab = "y", main = "gevcdn.bootstrap")
matplot(x, cbind(t(apply(CI$location.bootstrap, 1, quantile,
                              p = c(0.025, 0.975)))), type = c("l", "b", "b"), pch = 20,
        lwd = 3, col = c("black", rep("green", 2)), xlab = "x",
        ylab = "location", main = "location CI")
matplot(x, cbind(t(apply(CI$scale.bootstrap, 1, quantile,
                              p = c(0.025, 0.975)))), type = c("l", "b", "b"), pch = 20,
        lwd = 3, col = c("black", rep("orange", 2)), xlab = "x",
        ylab = "scale", main = "scale CI")
matplot(x, cbind(t(apply(CI$shape.bootstrap, 1, quantile,
                              p = c(0.025, 0.975)))), type = c("l", "b", "b"), pch = 20,
        lwd = 3, col = c("black", rep("brown", 2)), xlab = "x",
        ylab = "shape", main = "shape CI")


df.boot <- data.frame(year = x , obs = y,
                     loc = t(apply(CI$location.bootstrap, 1, quantile,
                                p = c(0.025, 0.975))),
                     scale = t(apply(CI$scale.bootstrap, 1, quantile,
                                  p = c(0.025, 0.975))),
                  t(apply(CI$shape.bootstrap, 1, quantile,
                          p = c(0.025, 0.975))))
gg.boot <- ggplot(df.boot, aes(x = year, y = obs)) +
  geom_line() + geom_point() +
  geom_line(aes(y = obs), col = "red", linetype = "dashed") +
  geom_line(aes(y = loc.2.5.), col = "blue") +
  geom_line(aes(y = loc.97.5.), col = "green") +
  labs(title = "GEV-CDN bagging with early stopping and 2 hidden nodes",
       y = "") +
  theme_piss()
gg.boot

