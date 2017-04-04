
functions {
real gev_log(real y, int t, real mu0, real mu1, real logsig, real xi){
    real mu; 
    real z;
    mu = mu0 + mu1 * t;
    z = 1 + (y - mu) * xi / exp(logsig);
  return -logsig - (1 + 1/xi)*log(z) - pow(z,-1/xi);
 }
}

data {
  int<lower=0> n;
  vector[n] y;
  int t;
}
parameters {
 real xi;
 real logsig;
 real mu0;
 real mu1;
}

// transformed parameters{
//   vector[num_elements(y)] mu;
//   for(j in 1:num_elements(y) ){
//     mu[j] = mu0 + mu1 * t[j];
//   }
// }
model {
 xi ~ normal(-0.1, 5);  # second arg is SD !!!!!! 
  logsig ~ normal(0, 5);
   mu0 ~ normal(30, 5);
   mu1 ~ normal(0, 5);
  
    y ~ gev(mu0, mu1, logsig, xi);
   
}
