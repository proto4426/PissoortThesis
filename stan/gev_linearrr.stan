functions {

  real gev_log (vector y, real mu, real sigma, real xi) { 
    vector[num_elements(y)] z; 
    vector[num_elements(y) + 1] lp; 
    real inv_xi; 
    real inv_xi_p1; 
    real neg_inv_xi; 
    z = ((1 + y) - mu) * (xi / sigma); 
    inv_xi = inv(xi); 
    inv_xi_p1 = 1 + inv_xi; 
    neg_inv_xi = -inv_xi; 
    for (n in 1:num_elements(y)) 
      lp[n] =  inv_xi_p1 * log(z[n]) + pow(z[n],neg_inv_xi); 
    lp[num_elements(y) + 1] = num_elements(y) * log(sigma); 
    return -sum(lp); 
  } 

}

data {
  int<lower=0> n;
  real y[n];
}

transformed data {
  real min_y;
  real max_y;
  real sd_y;
  min_y = min(y);
  max_y = max(y);
  sd_y = sd(y);
}
parameters {
  real<lower=-0.5,upper=0.5> xi;
  real<lower=0> sigma;
  // location has upper/lower bounds depending on the value of xi
  real<lower=if_else( xi < 0, min_y + sigma / xi, negative_infinity() ),
       upper=if_else( xi > 0, positive_infinity(), max_y + sigma / xi )> mu;
}
model {
  # priors
  sigma ~ normal(sd_y,100);
  #mu ~ uniform
  #xi ~ uniform(-.5,.5);
  xi ~ uniform(-.5,.5);

  for(i in 1:N) {
    y[i] ~ gev(mu, sigma, xi);  
  }
}


