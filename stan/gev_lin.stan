functions {
  real gev_lpdf(real y, real tt, real mu0, real mu1, real sigma, real xi) {
    real ystar;
    real xy_p1;
    real mu;
    mu = mu0 + mu1 * tt;
    ystar = (y - mu) / sigma;
    xy_p1 = xi * ystar + 1;
    return -log(sigma) - (1 + 1 / xi) * log(xy_p1) - (xy_p1) ^ (-1 / xi);
  }
}

data {
  int n;
  vector[n] y;
  vector[n] tt;
}

transformed data {
  real sd_y;
  real mean_y;
  mean_y = mean(y);
  sd_y = sd(y);

}

parameters {
  real mu0;
  real mu1;
  real<lower = 0> sigma;
  real xi;
}


model {
   sigma ~ normal(sd_y, 10);
  mu0 ~ normal(30, 10);
  xi ~ normal(0, 10);
  mu1 ~ normal(0, 10);

  for (i in 1:n) {
    y[i] ~ gev(tt[i], mu0, mu1, sigma, xi);
  }
}
