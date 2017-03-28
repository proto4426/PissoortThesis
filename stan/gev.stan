functions {
  real gev_lpdf(real y, real mu, real logsig, real xi) {
    real ystar;
    real xy_p1;

    ystar = (y - mu) / exp(logsig);
    xy_p1 = xi * ystar + 1;
    return -logsig - (1 + 1 / xi) * log(xy_p1) - (xy_p1) ^ (-1 / xi);
  }
}

data {
  int n;
  vector[n] y;
}


transformed data {
  real sd_y;
  real mean_y;
  mean_y = mean(y);
  sd_y = sd(y);
}

parameters {
  real mu;
  real logsig;
  real xi;
}


model {
   logsig ~ normal(log(sd_y), 10);
  mu ~ normal(mean(y), 10);
  xi ~ normal(0, 10);

  for (i in 1:n) {
    y[i] ~ gev(mu, logsig, xi);
  }
}



