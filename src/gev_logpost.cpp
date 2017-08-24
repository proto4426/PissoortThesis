// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// using namespace Rcpp;
// using namespace arma;

#include <Rcpp.h>
using namespace Rcpp;

// Miscellaneous functions

//' @export
// [[Rcpp::export]]
bool any_nonpos(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(x <= 0)) ;
}


// // [[Rcpp::export]]
// double cpp_gev_loglik0(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
//   if (x[1] <= 0)
//     return R_NegInf ;
//   NumericVector data = ss["data"] ;
//   NumericVector sdat = (data - x[0]) / x[1] ;
//   NumericVector zz = 1 + x[2] * sdat ;
//   if (any_nonpos(zz))
//     return R_NegInf ;
//   int m = ss["m"] ;
//   double val = -m * log(x[1]) ;
//   if (std::abs(x[2]) > 1e-6) {
//     val = val - (1 + 1 / x[2]) * sum(log(zz)) - sum(pow(zz, (-1 / x[2]))) ;
//   } else {
//     double sum_gev = ss["sum_gev"] ;
//     double t1, t2, sdatj, temp ;
//     double t0 = (sum_gev - m * x[0]) / x[1] ;
//     double tot = 0.0;
//     double tot2 = 0.0;
//     for(int j = 0; j < m; ++j) {
//       sdatj = sdat[j] ;
//       temp = 0.0 ;
//       for(int i = 1; i < 5; ++i) {
//         t1 = pow(sdatj, i) ;
//         t2 = (i * sdatj - i - 1) ;
//         tot += pow(-1.0, i) * t1 * t2 * pow(x[2], i) / i / (i + 1) ;
//         temp += pow(-1.0, i) * pow(sdatj, i + 1) * pow(x[2], i) / (i + 1) ;
//       }
//       tot2 += exp(-sdatj - temp) ;
//     }
//     val = val - t0 - tot - tot2 ;
//   }
//   return val ;
// }

// The 2ND PARAMETER IS LOGSIGMA !!!!
//' @export
// [[Rcpp::export]]
double cpp_gev_loglik(const NumericVector x, const NumericVector data) {
  if (x[1] <= 0)
    return R_NegInf ;
  // NumericVector data = ss["data"] ;
  NumericVector sdat = (data - x[0]) / log(x[1]) ;
  NumericVector zz = 1 + x[2] * sdat ;
  // Rcout << "zz are " << zz << std::endl ;
  if (any_nonpos(zz))
     return R_NegInf ;
  int m = data.size() ;
  double val = -m * log(log(x[1])) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - sum(pow(zz, (-1 / x[2]))) ;
  }
  else { }
}


// // For a user-definer prior
// // [[Rcpp::export]]
// double gev_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
//   SEXP prior_ptr = pars["prior"] ;
//   // Unwrap pointer to the log-prior function.
//   typedef double (*priorPtr)(const Rcpp::NumericVector& x,
//                   const Rcpp::List& ppars) ;
//   Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
//   priorPtr priorfun = *xpfun ;
//   return cpp_gev_loglik(x, pars) + priorfun(x, pars) ;
// }
//' @export
// [[Rcpp::export]]
double gev_logpost(const NumericVector& x,
                   const NumericVector data){//, const NumericVector priorsd ) {
  // Unwrap pointer to the log-prior function.

  // double mu = x[0] ;
  // double logsig = x[1] ;
  // double xi = x[2] ;
  //
  // double lprior ;
  // lprior = log(dnorm(1, mu, priorsd[0])) ;
  // lprior = lprior + log(dnorm(1, logsig, priorsd[1])) ;
  // lprior = lprior + log(dnorm(1, xi, priorsd[2])) ;
  return cpp_gev_loglik(x, data) ; //+ lprior ;
}

//' @export
// [[Rcpp::export]]
List gibbs_statioCpp(NumericVector start, NumericVector data,
                     int iter, NumericVector propsd,
                     bool verbose ) {

  // Initialize storage
  NumericMatrix  out(iter, start.size() ) ;
  if(verbose == TRUE) Rcout << "out init " << out << std::endl ;

  // Starting values
  out(0, _) = start ;

  if(verbose == TRUE) Rcout << "out starts= " << out << std::endl;
  if(verbose == TRUE) Rcout << "starts= " << start << std::endl;

  double lpost_old = gev_logpost(start, data) ;
  if(verbose == TRUE)  Rcout << "lpost_old= " << lpost_old << std::endl;

  // Acceptances rates
  NumericMatrix acc_rates(iter, propsd.size()) ;
  NumericVector mean_acc_rates(propsd.size()) ;

  // Gibbs sampler -- big loop
  for (int t = 0; t < iter; t++) {


    // For mu
    NumericVector prop1 = Rcpp::rnorm(1, out(t, 0), propsd[0] ) ;
    if(verbose == TRUE)  Rcout << "prop1= " << prop1[0] << std::endl;


    double new1par1 = prop1[0] ;
    double new1par2 = out(t,1);
    double new1par3 = out(t,2) ;
    NumericVector new_param1 =  NumericVector::create(new1par1, new1par2, new1par3);
    if(verbose == TRUE)  Rcout << " new_param1 init= " << new_param1 << std::endl;
    // newparam1 = NumericVector::create(prop1[0], out[t,1], out[t,2]);
    // Rcout << " new_param1 then= " << new_param1 << std::endl;

    double lpost_prop1 = gev_logpost(new_param1, data) ;
    if(verbose == TRUE)  Rcout << " lpost_prop1= " << lpost_prop1 << std::endl;

    double r1 = exp(lpost_prop1 - lpost_old) ;
    if(verbose == TRUE) Rcout << "r1= " << r1 << std::endl;

    NumericVector rcomp1 = runif(1) ;
    if(verbose == TRUE)  Rcout << "rcomp1= " << rcomp1 << std::endl;
     if(r1 >  rcomp1[0]){
       out(t+1,0) = prop1[0] ;
       lpost_old = lpost_prop1 ;
     }
     else {
       out(t+1,0) = out(t,0) ;
     }
     acc_rates(t,0) = std::min(r1,1.) ;


     // For logsig
     NumericVector prop2 = rnorm(1, out(t, 1), propsd[1]) ;
     if(verbose == TRUE)  Rcout << "prop2= " << prop2[0] << std::endl;

     double new2par1 = out(t+1,0) ;
     double new2par2 = prop2[0];
     double new2par3 = out(t,2) ;
     NumericVector new_param2 = NumericVector::create(new2par1, new2par2, new2par3);
     if(verbose == TRUE) Rcout << " new_param2 init= " << new_param2 << std::endl;

     double lpost_prop2 = gev_logpost(new_param2, data) ;

     double r2 = exp(lpost_prop2 - lpost_old) ;
     if(verbose == TRUE) Rcout << "r2= " << r2 << std::endl;

     NumericVector rcomp2 = runif(1) ;
     if(r2 > rcomp2[0]){
       out(t+1,1) = prop2[0] ;
       lpost_old = lpost_prop2 ;
     }
     else {
       out(t+1,1) = out(t,1) ;
     }
     acc_rates(t,1) = std::min(r2, 1.) ;


     // For xi
     NumericVector prop3 = rnorm(1, out(t, 2), propsd[2]) ;
     if(verbose == TRUE)  Rcout << "prop3= " << prop3[0] << std::endl;

     double new3par1 = out(t+1,0) ;
     double new3par2 = out(t+1,1);
     double new3par3 = prop3[0] ;
     NumericVector new_param3 = NumericVector::create(new3par1, new3par2, new3par3) ;
     if(verbose == TRUE) Rcout << " new_param3 init= " << new_param3 << std::endl;

     double lpost_prop3 = gev_logpost(new_param3, data) ;

     double r3 = exp(lpost_prop3 - lpost_old) ;
     if(verbose == TRUE) Rcout << "r3= " << r3 << std::endl;

     NumericVector rcomp3 = runif(1) ;
      if(r3 > rcomp3[0]){
       out(t+1,2) = prop3[0] ;
       lpost_old = lpost_prop3 ;
     }
     else {
       out(t+1,2) = out(t,2) ;
     }
     acc_rates(t,2) = std::min(r3, 1.) ;


     if(verbose == TRUE)  Rcout << "lpost_old end loop= " << lpost_old << std::endl;

     mean_acc_rates += acc_rates(t, _) ;
  }

  // MEan acceptances rates
  mean_acc_rates = mean_acc_rates / iter ;


  return Rcpp::List::create(Rcpp::Named("out.chain")=out,
                            Rcpp::Named("acc.rates") = acc_rates,
                            Rcpp::Named("mean.acc.rates") = mean_acc_rates) ;
                           // Rcpp::Named("iter")=iter);
}
