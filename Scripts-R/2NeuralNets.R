setwd('/home/piss/Documents/Extreme/R resources/IRM')
load("data1.Rdata")

library(PissoortThesis)

library(GEVcdn)

# fit <- gevcdn.fit(as.matrix(seq(1, length(max_years$data))), as.matrix(max_years$data))
# gevcdn.evaluate(as.matrix(seq(1, length(max_years$data))), fit)
#
# fit <- gevcdn.fit(as.matrix(seq(1, length(max_all))), as.matrix(max_all))
# ff <- gevcdn.evaluate(as.matrix(seq(1, length(max_all))), fit)



## 1) Define the hierarchy of models of increasing complexity

models <- list()
# Stationary model
models[[1]] <- list(Th = gevcdn.identity,
                    fixed = c("location", "scale", "shape"))
# Linear models
models[[2]] <- list(Th = gevcdn.identity, fixed = c("shape","scale"))
models[[3]] <- list(Th = gevcdn.identity, fixed = c("shape"))
models[[4]] <- list(Th = gevcdn.identity)

# Nonlinear, 1 or 2 hidden nodes
models[[5]] <- list(n.hidden = 1, Th = gevcdn.logistic, fixed = c("shape", "scale"))
models[[6]] <- list(n.hidden = 2, Th = gevcdn.logistic, fixed = c("shape", "scale"))
models[[7]] <- list(n.hidden = 1, Th = gevcdn.logistic, fixed = c("shape"))
models[[8]] <- list(n.hidden = 2, Th = gevcdn.logistic, fixed = c("shape"))
models[[9]] <- list(n.hidden = 1, Th = gevcdn.logistic)
models[[10]] <- list(n.hidden = 2, Th = gevcdn.logistic)





x <- as.matrix(seq(1, length(max_years$data)))
y <- as.matrix(max_years$data)
## 2) Fit the models and retrieve the weights

weights.models <- list()
for(i in seq_along(models)){
  weights.models[[i]] <- gevcdn.fit(x = x, y = y, n.trials = 1,
                                    n.hidden = models[[i]]$n.hidden,
                                    Th = models[[i]]$Th,
                                    fixed = models[[i]]$fixed)
}
# Printed outputs correspond to the value of the optimized ?gevcdn.cost (see also ?optim())
#for each model in this loop.



## Select best model

models.AICc <- round(sapply(weights.models, attr, which = "AICc"), 3)
# Comparing the AICc, we confirm that shape parameter must be held fixed.
# But last model seems also good...
models.BIC <- round(sapply(weights.models, attr, which = "BIC"), 3)
# Clear evidence for the 5th model (simple linear trend in location parameter)
# BIC penalizes more the more complex models --> pasimony
weights.best <- weights.models[[which.min(models.AICc)]]
parms.best <- gevcdn.evaluate(x, weights.best)


#### Compare nested models : Deviance statistics

## stationary vs linear trend
nll1 <- attr(weights.models[[1]], "NLL")
nll5 <- attr(weights.models[[5]], "NLL")
pchisq( 2 *( (-nll5) - (-nll1) ), lower.tail = F,
        df = attr(weights.models[[5]], "k") - attr(weights.models[[1]], "k"))
# ~exactly same result as previously done : linear trend clearly significant

## we can easily guess the output produced for the other tests

q.best <- sapply(c(0.1, 0.5, 0.9), qgev,
                 location = parms.best[,"location"],
                 scale = parms.best[,"scale"],
                 shape = parms.best[,"shape"])


############ Bagging ################

n.boot <- 10 # Number of boostrapped iterations

## Do it through parallel computing
weights.on <- gevcdn.bag(x = x, y = y, iter.max = 100,
                         iter.step = 10, n.bootstrap = n.boot,
                         n.hidden = 2)
parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)

# parms <- list()
# for (i in 1:n.boot){
#   parms[[i]] <- apply(parms.on[[i]],2, mean)
# }
parms <- apply(parms.on, 2, FUN = mean)  ## Try with this

parms <- matrix(0, nrow = nrow(parms.on[[1]]), ncol = 3)
for (i in 1:n.boot){
  parms <- parms + as.matrix(parms.on[[i]])
}
parms <- parms / length(parms.on)


q <- t(sapply(max_years$data, quantile, probs = c(.1, .5, .9)))

q.10.on <- q.50.on <- q.90.on <- c()
for(i in seq_along(parms.on)){
  q.10.on <- cbind(q.10.on, VGAM::qgev(p = 0.1,
                                 location = parms.on[[i]][,"location"],
                                 scale = parms.on[[i]][,"scale"],
                                 shape = parms.on[[i]][,"shape"]))
  q.50.on <- cbind(q.50.on, VGAM::qgev(p = 0.5,
                                 location = parms.on[[i]][,"location"],
                                 scale = parms.on[[i]][,"scale"],
                                 shape = parms.on[[i]][,"shape"]))
  q.90.on <- cbind(q.90.on, VGAM::qgev(p = 0.9,
                                 location = parms.on[[i]][,"location"],
                                 scale = parms.on[[i]][,"scale"],
                                 shape = parms.on[[i]][,"shape"]))
}
## Plot data and quantiles
matplot(cbind(y, q, rowMeans(q.10.on), rowMeans(q.50.on),
              rowMeans(q.90.on)), type = c("b", rep("l", 6)),
        lty = c(1, rep(c(1, 2, 1), 2)),
        lwd = c(1, rep(c(3, 2, 3), 2)),
        col = c("red", rep("orange", 3), rep("blue", 3)),
        pch = 19, xlab = "x", ylab = "y",
        main = "gevcdn.bag (early stopping on)")




# More relevant to do this with rainfall data !?
#############
datap50 <- read.csv('P50_Uccle_1880.csv',sep="")
#################



