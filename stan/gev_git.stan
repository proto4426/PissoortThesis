// // functions {
// //   real gev_lpdf(real y, real mu, real logsig, real xi) {
// //     real ystar;
// //     real xy_p1;
// //
// //     ystar = (y - mu) / exp(logsig);
// //     xy_p1 = xi * ystar + 1;
// //     return -logsig - (1 + 1 / xi) * log(xy_p1) - (xy_p1) ^ (-1 / xi);
// //   }
// // }
// //
// // data {
// //   int n;
// //   vector[n] y;
// // }
// //
// // transformed data {
// //   real mean_y;          real sd_y;
// //   real min_y;         real max_y;
// //
// //   mean_y = mean(y);     sd_y = sd(y);
// //   min_y = min(y);     max_y = max(y);
// // }
// //
// // parameters {
// //   real logsig;
// //   real xi;
// //   real <lower=if_else( xi < 0, min_y + exp(logsig)/xi, negative_infinity() ),
// //     upper=if_else( xi > 0, positive_infinity(), max_y + exp(logsig)/xi )>
// //     mu;
// // }
// //
// // model {
// //   xi ~ normal(0, 10);
// //
// //   for (i in 1:n) {
// //     y[i] ~ gev(mu, logsig, xi);
// //   }
// // }
//
//
// functions {
//
//  real gev_lpdf (vector y, real mu, real logsig, real xi) {
//   vector[rows(y)] z;
//   vector[rows(y) + 1] lp;
//   real inv_xi;
//   real inv_xi_p1;
//   real neg_inv_xi;
//   z = ((1 + y) - mu) * (xi / exp(logsig) );
//   inv_xi = inv(xi);
//   inv_xi_p1 = 1 + inv_xi;
//   neg_inv_xi = -inv_xi;
//   for (n in 1:rows(y))
//     lp[n] =  inv_xi_p1 * log(z[n]) + pow(z[n],neg_inv_xi);
//   lp[rows(y) + 1] = rows(y) * logsig;
//   return -sum(lp);
// }
// }
//
// data {
//   int<lower=0> n;
//   vector[n] y;
// }
//
// transformed data {
//   real min_y;
//   real max_y;
//   min_y = min(y);
//   max_y = max(y);
// }
// parameters {
//   real<lower=-.5, upper=0> xi;
//   real logsig;
//   // location has upper/lower bounds depending on the value of xi
//   real mu;
// }
//
// model {
//    xi ~ normal(-0.1, 5);  # second arg is SD !!!!!!
//   logsig ~ normal(0, 5);
//    mu ~ normal(30, 5);
//    y ~ gev(mu,logsig,xi);
// }


functions {

  real gev_log (real y, real mu, real sigma, real xi){
    real z;
   # real<lower=-0.5,upper=0.15> xi;
    z = 1 + (y - mu) * xi / sigma;
    return -log(sigma) - (1 + 1/xi)*log(z) - pow(z,-1/xi);
  }

}

data {
  int<lower=0> n;
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
  // location has upper/lower bounds depending on the value of xi
  real<lower=if_else( xi < 0, min_y + sigma / xi, negative_infinity() ),
       upper=if_else( xi > 0, positive_infinity(), max_y + sigma / xi )> mu;
}
model {
  # priors on component parameters
  sigma ~ gamma(.0001,.0001);
  #mu ~ uniform
  xi ~ uniform(-.5,.15);

  for(i in 1:n) {
    y[i] ~ gev(mu, sigma, xi);
  }
}
