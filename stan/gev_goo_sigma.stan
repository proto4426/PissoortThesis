data {
  int<lower=0> n;
  vector[n] t;
  real y[n];
}

transformed data {
  real min_y;
  real max_y;
  min_y <- min(y);
  max_y <- max(y);
}
parameters {
  real xi;
  real<lower=0> sigma;
  real<lower=if_else( xi > 0, min_y + sigma / xi, negative_infinity() ),
       upper=if_else( xi > 0, positive_infinity(), max_y + sigma / xi )> mu;
}
model {
  real z;
  for(i in 1:n) {
    z = 1 + (y[i] - mu) * xi / sigma;
    increment_log_prob(-log(sigma) - (1 + 1/xi)*log(z) - pow(z,-1/xi));
  }
  // priors on xi, sigma, and mu

  // sigma ~ normal(log(sd_y), 50);
  // #mu ~ normal(mean_y, 50);
  // #xi ~ normal(0, 50);
  // xi ~ uniform(-1, 0.5);
}
