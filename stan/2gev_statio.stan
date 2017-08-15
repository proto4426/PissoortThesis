functions {
real gev_log (vector y, real mu, real logsig, real xi){
   vector[num_elements(y)] z;
   vector[num_elements(y) ] res;
   // real<lower=0> z[n];
   # real<lower=-0.5,upper=0.15> xi;
   z = ((1 + y) - mu) * (xi / exp(logsig));
   for (i in 1:num_elements(y))
      z[i] = if_else(z[i]<=0, 0.0000001, z[i] );
   for(j in 1:num_elements(y) )
     res[j] = -logsig - (1 + 1/xi)*log(z[j]) - pow(z[j], -1/xi);
  return -sum(res);
  }
}
data {
  int<lower=0> n;
  vector[n] y;
  #int y[n];
}
transformed data {
  real min_y;
  real max_y;
  min_y = min(y);
  max_y = max(y);
}
parameters {
 real xi;
 real logsig;
  // location has upper/lower bounds depending on the value of xi
 real<lower=if_else( xi < 0, min_y + logsig / xi, negative_infinity() ),
      upper=if_else( xi > 0, positive_infinity(), max_y + logsig / xi )>
      mu;
       // vector[3] theta;
       // real theta = [mu, sigma, xi ];
 #vector[3] theta;
}
model {
 #xi ~ normal(-0.1, 5);  # second arg is SD !!!!!!
 # sigma ~ normal(1, 10);
 # priors on component parameters
  logsig ~ normal(0, 5);
   mu ~ normal(30, 10);
  xi ~ uniform(-1,0);
  # theta += normal_lpdf(y | real[30, 0, 0], real[10, 10, 10])

    y ~ gev(mu, logsig, xi);

  #y ~ gev(mu, sigma, xi);
   // target += normal_lpdf(y | 1, 10);
   // target += normal_lpdf(y | 30, 10);
   // target += normal_lpdf(y | -0.1, 10);
}


