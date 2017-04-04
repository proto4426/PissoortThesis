functions {

 real gev_lpdf (vector y, real mu, real logsig, real xi) {
  vector[rows(y)] z;
  vector[rows(y) + 1] lp;
  real inv_xi;
  real inv_xi_p1;
  real neg_inv_xi;
  z = ((1 + y) - mu) * (xi / exp(logsig) );
  inv_xi = inv(xi);
  inv_xi_p1 = 1 + inv_xi;
  neg_inv_xi = -inv_xi;
  for (n in 1:rows(y))
    lp[n] =  inv_xi_p1 * log(z[n]) + pow(z[n],neg_inv_xi);
  lp[rows(y) + 1] = rows(y) * logsig;
  return -sum(lp);
} 
}

data {
  int<lower=0> N;
  vector[N] y;
}

transformed data {
  real min_y;
  real max_y;
  min_y = min(y);
  max_y = max(y);
}
parameters {
  real<lower=-.5, upper=0> xi;
  real logsig;
  // location has upper/lower bounds depending on the value of xi
  real mu;
}

model {
   xi ~ normal(-0.1, 5);  # second arg is SD !!!!!! 
  logsig ~ normal(0, 5);
   mu ~ normal(30, 5);
   y ~ gev(mu,logsig,xi); 
}

