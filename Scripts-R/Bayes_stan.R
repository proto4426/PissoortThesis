#setwd('/home/piss/Documents/Extreme/R resources/IRM')
setwd('/home/piss/PissoortRepo/PissoortThesis/stan')

library("rstantools")

load("/home/piss/Documents/Extreme/R resources/IRM/data1.Rdata")

options(mc.cores=parallel::detectCores()) # all available cores
# can be used without needing to manually specify the cores argument.

library("rstan")      
library(bayesplot)
library(mvtnorm)

library(PissoortThesis)


#######
fn <- function(par, data) -log_post0(par[1], par[2], par[3], data)
param <- c(mean(max_years$df$Max),log(sd(max_years$df$Max)), 0.1 )
opt <- nlm(fn, param, data = max_years$data,
           hessian=T, iterlim = 1e5) 

start0 <- list() ;  k <- 1
while(k < 5) { # starting value is randomly selected from a distribution
  # that is overdispersed relative to the target
  sv <- as.numeric(rmvnorm(1, opt$estimate, 50 * solve(opt$hessian)))
  svlp <- log_post0(sv[1], sv[2], sv[3], max_years$data)
  if(is.finite(svlp)) {
    start0[[k]] <- as.list(sv) ;  names(start0[[k]]) <- c("mu", "logsig","xi") 
    k <- k + 1
  } }
######

fit_stan <- stan(file = 'gev.stan', data = list(n = length(max_years$data),
                                                y = max_years$data), 
                 iter = 2000, chains = 4, warmup = 0,# init = rev(start0), 
                 cores = 8, verbose = T, control = list(adapt_delta = .9))
fit_stan
summary(fit_stan)
pairs(fit_stan)
sampler_par <- get_sampler_params(fit_stan, inc_warmup = TRUE)
summary(do.call(rbind, sampler_par), digits = 2)
lapply(sampler_par, summary, digits = 2)

lookup(Inf)



tt <- ( min(max_years$df$Year):max(max_years$df$Year) -
          mean(max_years$df$Year) ) / length(max_years$data)
start <- list() ; k <- 1
while(k < 5) { 
  sv <- as.numeric(rmvnorm(1, opt$par, 50 * solve(opt$hessian)))
  svlp <- log_post1(sv[1], sv[2], sv[3], sv[4], max_years$data)
  print(svlp)
  if(is.finite(svlp)) {
    start[[k]] <- sv;  names(start[[k]]) <- c("mu0", "mu1", "sigma","xi") 
    k <- k + 1
  }
}

fit_lin <- stan(file = 'gev_lin.stan', data = list(n = length(max_years$data),
                                                  y = max_years$data,
                                                  tt = tt), 
              init=start, iter = 2000, chains = 4, warmup = 100, cores = 8,
               control = list(adapt_delta = .99))
fit_lin
sampler_params <- get_sampler_params(fit_lin, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
lapply(sampler_params, summary, digits = 2)

arrayfit <- as.array(fit_stan) 
mcmc_trace(arrayfit, pars = c("mu", "logsig", "xi"),
           facet_args = list(ncol = 1, strip.position = "left"))

test_fun <-' functions {
  real gev_log (vector y, real mu, real logsig, real xi) { 
    vector[num_elements(y)] z; 
    vector[num_elements(y) + 1] lp; 
    real inv_xi; 
    real inv_xi_p1; 
    real neg_inv_xi; 
    z = ((1 + y) - mu) * (xi / exp(logsig)); 
    inv_xi = inv(xi); 
    inv_xi_p1 = 1 + inv_xi; 
    neg_inv_xi = -inv_xi; 
    for (n in 1:num_elements(y)) 
      lp[n] =  inv_xi_p1 * log(z[n]) + pow(z[n],neg_inv_xi); 
    lp[num_elements(y) + 1] = num_elements(y) * logsig; 
    return -sum(lp); 
  } }
  model {}
'
test_fun <- ' functions { 
  real gg( int y) {
   return 1e15+y; 
  } 
}   
  model {}
'
expose_stan_functions(stanc(model_code = test_fun))
gev_log(max_years$data, 30, 0.7, -.22)
gg(20)
sm <- stan_model(model_code = test_fun)
optimizing(sm, data = list(y = max_years$data,
                           N = length(max_years$data)), hessian = T)


############# Diagnostics + visualization with bayesplot package 

library(bayesplot) # Mainly for use with Stan 


mcmc_intervals(param.chain, pars = c("logsig", "xi"))
# Bring back mu
mcmc_areas(
  param.chain, 
  pars = c("mu1", "logsig", "xi"),
  prob = 0.9, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)
color_scheme_set("brightblue")
grid.arrange( mcmc_dens(param.chain, pars = c("mu", "mu1")),
              mcmc_dens(param.chain, pars = c("logsig", "xi")),
              nrow = 2)
str(gibbs.trend$out.ind)
str(param.chain)
un.out <- unlist(gibbs.trend$out.ind,use.names = T, recursive = F)
matrix(unlist(gibbs.trend$out.ind), ncol = 4, byrow = F)
unlst <- do.call(rbind, gibbs.trend$out.ind )

array.post <- array(unlist(t(gibbs.trend$out.ind)), dim = c(4000/4+1, 4, 3 ))
# dimnames = list(c(NULL, 
#           c("mu", "mu1", "logsig", "xi"),
#           c("chain:1", "chain:2", "chain:3", "chain:4"))))
dim(array.post)
dimnames(array.post) <- list(iterations = NULL, 
                             parameters = c("mu", "mu1", "logsig", "xi"),
                             chains = c("chain:1", "chain:2", "chain:3"))

array.post <- aperm(array.post, perm = c("iterations", "chains", "parameters"))
dimnames(array.post)
str(array.post)
color_scheme_set("green")
mcmc_hist_by_chain( array.post, 
                    pars = c("mu", "mu1", "logsig", "xi")) 


mcmc_dens_overlay(array.post, pars = c("mu", "mu1"))


color_scheme_set("mix-blue-red")
mcmc_trace(posterior, pars = c("wt", "sigma"), 
           facet_args = list(ncol = 1, strip.position = "left"))
