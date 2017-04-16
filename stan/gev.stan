functions {
  real gev_lpdf(real y, real mu, real logsig, real xi) {
    real ystar;
    real xy_p1;

    ystar = (y - mu) / exp(logsig);
    xy_p1 = xi * ystar + 1;
    // if (xy_p1 <= 0) return 1e8;
    // else return -logsig - (1 + 1 / xi) * log(xy_p1) - (xy_p1) ^ (-1 / xi);
    return -logsig - (1 + 1 / xi) * log(xy_p1) - (xy_p1) ^ (-1 / xi);
  }
}

data {
  int<lower=0> n;
  vector[n] y;
 }

transformed data {
  real mean_y;          real sd_y;
  real min_y;         real max_y;
  
  mean_y = mean(y);     sd_y = sd(y);
  min_y = min(y);     max_y = max(y);
}

parameters {
  real logsig;
  real xi;
  real
  <lower=if_else( xi < 0, min_y + exp(logsig) / xi, negative_infinity() ),
    upper=if_else( xi > 0, positive_infinity(), max_y + exp(logsig) / xi )>
       mu;
}


model {
  logsig ~ normal(log(sd_y), 50);
 # mu ~ normal(mean_y, 50);
  xi ~ normal(0, 50);

   # y ~ gev(mu, logsig, xi);  

  for (i in 1:n) {
    y[i] ~ gev_lpdf(mu, logsig, xi);
  }
}



