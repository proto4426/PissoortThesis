#include <Rcpp.h>
using namespace Rcpp;


//' @export
// [[Rcpp::export]]
bool any_nonpos(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(x <= 0)) ;
}

// The 2ND PARAMETER IS LOGSIGMA !!!!
//' @export
// [[Rcpp::export]]
double cpp_gev_loglik(const NumericVector x,
                      const NumericVector data) {
  NumericVector sdat = (data - x[0]) / exp(x[1]) ;
  NumericVector zz = 1 + x[2] * sdat ;
 // Rcout << "zz are " << zz << std::endl ;
  if (any_nonpos(zz) || x[1]<=1e-10)
     return R_NegInf ;
  int m = data.size() ;
  double val = -m * x[1] ;
 // Rcout << "val is " << val << std::endl;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - sum(pow(zz, (-1 / x[2])));
  }
  else { // Gumbel case
    val = val - sum(sdat) - sum(exp(-sdat)) ;
    }
  return val ;
}

// LOGsigma
// [[Rcpp::export]]
double cpp_gev_loglik_other(const NumericVector x,
                            const NumericVector data){
  // if (x[1] <= 0)
  //   return R_NegInf ;
  NumericVector z = 1 + (data - x[0]) * (x[2] / exp(x[1]) );
  int N = data.size() ;
  // Rcout << "THE VECTOR z is " << z << std::endl ;
  NumericVector lp(N+1);
  for (int n = 0; n < N; n++){
     if(z[n] < 0) z[n] = 0.000001 ;
     lp[n] =  (1+pow(x[2],-1)) * log(z[n]) + pow(z[n],-pow(x[2],-1));
    // Rcout << "Lp is " << lp[n] << std::endl;
  }
  lp[N + 1] = N * x[1];
  return -sum(lp);
}


//========= Stationary Model ==============================================


//' @export
// [[Rcpp::export]]
double gev_logpost(const NumericVector& x,
                   const NumericVector data){
  double mu = x[0] ;
  double logsig = x[1] ;
  double xi = x[2] ;

  double lprior ;  // Large Variance priors
  lprior = R::dnorm(mu, 0., 50., TRUE) ;
  lprior = lprior + R::dnorm(logsig, 0., 50., TRUE) ;
  lprior = lprior + R::dnorm(xi, 0., 10., TRUE) ;

  double loglik = cpp_gev_loglik(x, data) ;

  return loglik + lprior ;
}


//' @useDynLib PissoortThesis
//' @export
// [[Rcpp::export]]
NumericMatrix gibbs_statioCpp(NumericVector start, NumericVector data,
                     int iter, NumericVector propsd,
                     bool verbose ) {

  // Initialize storage
  NumericMatrix  out(iter+1, start.size() ) ;
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
    for (int t = 0; t<iter; t++){

    // For mu
    NumericVector prop1 = Rcpp::rnorm(1, out(t, 0), propsd[0] ) ;
    if(verbose == TRUE)  Rcout << "prop1= " << prop1[0] << std::endl;

    NumericVector new_param1 =  NumericVector::create(prop1[0] , out(t,1), out(t,2));
    if(verbose == TRUE)  Rcout << " new_param1 init= " << new_param1 << std::endl;

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


     // For LOGSIGMA
     NumericVector prop2 = rnorm(1, out(t, 1), propsd[1]) ;
     if(verbose == TRUE)  Rcout << "prop2= " << prop2[0] << std::endl;

     NumericVector new_param2 = NumericVector::create(out(t+1,0), prop2[0], out(t,2));
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

     NumericVector new_param3 = NumericVector::create(out(t+1,0), out(t+1,1), prop3[0]);
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
  Rcout << "mean_acc_rates FINAL= " << mean_acc_rates << std::endl;


  return out; //  Rcpp::NumericMatrix::create(

                            // Rcpp::Named("acc.rates") = acc_rates,
                            // Rcpp::Named("mean.acc.rates") = mean_acc_rates) ;
                           // Rcpp::Named("iter")=iter);
}



//========= Nonstationary Model ===========================================
// Linear trend on the location ===========
// =======================================================================


//' @export
// [[Rcpp::export]]
double gevNsta_loglik(const NumericVector x,
                      const NumericVector data,
                      NumericVector tt) {
  NumericVector mu = x[0] + x[1] * tt ;

  NumericVector sdat = (data - mu) / exp(x[2]) ;
  NumericVector zz = 1 + x[3] * sdat ;
  // Rcout << "zz are " << zz << std::endl ;
  if (any_nonpos(zz) || x[2]<=1e-10)
    return R_NegInf ;
  int m = data.size() ;
  double val = -m * x[2] ;
  // Rcout << "val is " << val << std::endl;
  if (std::abs(x[3]) > 1e-6) {
    val = val - (1 + 1 / x[3]) * sum(log(zz)) - sum(pow(zz, (-1 / x[3])));
  }
  else { // Gumbel case
    val = val - sum(sdat) - sum(exp(-sdat)) ;
  }
  return val ;
}
//' @export
// [[Rcpp::export]]
double gevNsta_lpost(const NumericVector& x,
                   const NumericVector data,
                   NumericVector tt,
                   NumericVector mnpr = NumericVector::create(30,0,0,0),
                   NumericVector sdpr = NumericVector::create(40,40,10,10)){
  double mu0 = x[0] ;
  double mu1 = x[1] ;
  double logsig = x[2] ;
  double xi = x[3] ;
  NumericVector theta = NumericVector::create(mu0, mu1, logsig, xi);

  NumericVector mu(tt.size()) ;
  mu = mu0 + mu1 * tt ;  // Linear model on mu

  double lprior ;  // Large Variance priors
  double lpriormu0 ;    double lpriormu1 ;
  double lpriorlogsig ;    double lpriorxi ;
  lpriormu0    = R::dnorm(mu0, mnpr[0], sdpr[0], TRUE) ;
  lpriormu1    = R::dnorm(mu1,mnpr[1], sdpr[1], TRUE) ;
  lpriorlogsig =  R::dnorm(logsig, mnpr[2], sdpr[2], TRUE) ;
  lpriorxi     =  R::dnorm(xi, mnpr[3], sdpr[3], TRUE) ;

  lprior = lpriormu0 + lpriormu1 + lpriorlogsig + lpriorxi ;

  double loglik  ;
  loglik = gevNsta_loglik(theta, data, tt);

  // for(int i = 0 ; i<tt.size(); i++){
  //   NumericVector new_x =  NumericVector::create(mu,  logsig, xi) ;
  //   gevNsta_loglik(new_x, data)
  // }

  return loglik + lprior ;
}


//' @useDynLib PissoortThesis
//' @export
// [[Rcpp::export]]
List gibbs_NstaCpp(const NumericVector start, const NumericVector data,
                            NumericVector tt,
                            int iter, NumericVector propsd,
                            NumericVector mnpr = NumericVector::create(30,0,0,0),
                            NumericVector sdpr = NumericVector::create(40,40,10,10),
                            bool verbose = TRUE) {

  // Initialize storage
  NumericMatrix  out(iter+1, start.size() ) ;
  if(verbose == TRUE) Rcout << "out init " << out << std::endl ;

  // Starting values
  out(0, _) = start ;

  if(verbose == TRUE) Rcout << "out starts= " << out << std::endl;
  if(verbose == TRUE) Rcout << "starts= " << start << std::endl;

  double lpost_old = gevNsta_lpost(start, data, tt, mnpr, sdpr) ;
  if(verbose == TRUE)  Rcout << "lpost_old= " << lpost_old << std::endl;

  // Acceptances rates
  NumericMatrix acc_rates(iter, propsd.size()) ;
  NumericVector mean_acc_rates(propsd.size()) ;


  // Gibbs sampler -- big loop
  for (int t = 0; t<iter; t++){

    // For mu0
    NumericVector prop1 = Rcpp::rnorm(1, out(t, 0), propsd[0] ) ;
    if(verbose == TRUE)  Rcout << "prop1= " << prop1[0] << std::endl;

    NumericVector new_param1 =  NumericVector::create(prop1[0] , out(t,1), out(t,2),out(t,3));
    if(verbose == TRUE)  Rcout << " new_param1 init= " << new_param1 << std::endl;

    double lpost_prop1 = gevNsta_lpost(new_param1, data, tt, mnpr, sdpr) ;
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


    // For mu1
    NumericVector prop2 = Rcpp::rnorm(1, out(t, 1), propsd[1] ) ;
    if(verbose == TRUE)  Rcout << "prop1= " << prop2[0] << std::endl;

    NumericVector new_param2 =  NumericVector::create(out(t+1,0), prop2[0], out(t,2),out(t,3));
    if(verbose == TRUE)  Rcout << " new_param1 init= " << new_param2 << std::endl;

    double lpost_prop2 = gevNsta_lpost(new_param2, data, tt, mnpr, sdpr) ;
    if(verbose == TRUE)  Rcout << " lpost_prop1= " << lpost_prop2 << std::endl;

    double r2 = exp(lpost_prop2 - lpost_old) ;
    if(verbose == TRUE) Rcout << "r1= " << r1 << std::endl;

    NumericVector rcomp2 = runif(1) ;
    if(verbose == TRUE)  Rcout << "rcomp1= " << rcomp1 << std::endl;
    if(r2 >  rcomp2[0]){
      out(t+1,1) = prop2[0] ;
      lpost_old = lpost_prop2 ;
    }
    else {
      out(t+1,1) = out(t,1) ;
    }
    acc_rates(t,1) = std::min(r2,1.) ;



    // For LOGSIGMA
    NumericVector prop3 = rnorm(1, out(t, 2), propsd[2]) ;
    if(verbose == TRUE)  Rcout << "prop2= " << prop3[0] << std::endl;

    NumericVector new_param3 = NumericVector::create(out(t+1,0), out(t+1,1), prop3[0], out(t,3));
    if(verbose == TRUE) Rcout << " new_param2 init= " << new_param3 << std::endl;

    double lpost_prop3 = gevNsta_lpost(new_param3, data, tt, mnpr, sdpr) ;

    double r3 = exp(lpost_prop3 - lpost_old) ;
    if(verbose == TRUE) Rcout << "r2= " << r2 << std::endl;

    NumericVector rcomp3 = runif(1) ;
    if(r3 > rcomp3[0]){
      out(t+1,2) = prop3[0] ;
      lpost_old = lpost_prop3 ;
    }
    else {
      out(t+1,2) = out(t,2) ;
    }
    acc_rates(t,2) = std::min(r3, 1.) ;


    // For xi
    NumericVector prop4 = rnorm(1, out(t, 3), propsd[3]) ;
    if(verbose == TRUE)  Rcout << "prop3= " << prop4[0] << std::endl;

    NumericVector new_param4 = NumericVector::create(out(t+1,0), out(t+1,1),out(t+1,2), prop4[0]);
    if(verbose == TRUE) Rcout << " new_param3 init= " << new_param4 << std::endl;

    double lpost_prop4 = gevNsta_lpost(new_param4, data, tt, mnpr, sdpr) ;

    double r4 = exp(lpost_prop4 - lpost_old) ;
    if(verbose == TRUE) Rcout << "r4= " << r3 << std::endl;

    NumericVector rcomp4 = runif(1) ;
    if(r4 > rcomp4[0]){
      out(t+1,3) = prop4[0] ;
      lpost_old = lpost_prop4 ;
    }
    else {
      out(t+1,3) = out(t,3) ;
    }
    acc_rates(t,3) = std::min(r4, 1.) ;


    if(verbose == TRUE)  Rcout << "lpost_old end loop= " << lpost_old << std::endl;

    mean_acc_rates += acc_rates(t, _) ;
  }

  // MEan acceptances rates
  mean_acc_rates = mean_acc_rates / iter ;
  // Rcout << "mean_acc_rates FINAL= " << mean_acc_rates << std::endl;

  // return out;
  return Rcpp::List::create( Rcpp::Named("out.ind") = out,
                             Rcpp::Named("mean.acc.rates") = mean_acc_rates
                               ) ;
}
