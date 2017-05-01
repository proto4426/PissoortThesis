#rsconnect::setAccountInfo(name, token, secret)  # FILL it ! 

library(shinystan)
my_sso <- launch_shinystan(fit_stan)




deploy_shinystan(my_sso, appName = "MyModel", account = "username")



