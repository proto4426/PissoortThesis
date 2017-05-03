#rsconnect::setAccountInfo(name, token, secret)  # FILL it !

# load(fit_stan.RData)

library(shinystan)
my_sso <- launch_shinystan(fit_stan) # Put you stan model in here.




deploy_shinystan(my_sso, appName = "MyModel", account = "username")



