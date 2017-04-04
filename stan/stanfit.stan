functions {

  real gev_log (real y, real mu, real sigma, real xi){
    real z;
   # real<lower=-0.5,upper=0.15> xi;
    z = 1 + (y - mu) * xi / sigma;
    return -log(sigma) - (1 + 1/xi)*log(z) - pow(z,-1/xi);  
  }

}

data {
  int<lower=0> N;
  real y[N];
}
transformed data {
  real min_y;
  real max_y;
  min_y = min(y);
  max_y = max(y);
}
parameters {
  real xi;
  real<lower=0> sigma;
  // location has upper/lower bounds depending on the value of xi
  real<lower=if_else( xi < 0, min_y + sigma / xi, negative_infinity() ),
       upper=if_else( xi > 0, positive_infinity(), max_y + sigma / xi )> mu;
}
model {
  # priors on component parameters
  sigma ~ gamma(.0001,.0001);
  #mu ~ uniform
  xi ~ uniform(-.5,.15);

  for(i in 1:N) {
    y[i] ~ gev(mu, sigma, xi);  
  }
}

