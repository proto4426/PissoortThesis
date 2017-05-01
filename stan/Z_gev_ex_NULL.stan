data {
  int<lower=0> n;
  #vector[n] t;
  real y[n];
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
  real<lower=if_else( xi < 0, min_y + sigma / xi, -1e15),
  upper=if_else( xi > 0, 1e15, max_y + sigma / xi )> mu;
}
model {
  real z[n];
  for(i in 1:n) {
    z[n] = 1 + (y[i] - mu) * xi / sigma;
    target += -log(sigma) - (1 + 1/xi)*log(z[n]) - pow(z[n],-1/xi);  
  }
  // priors on xi, sigma, and mu
}
