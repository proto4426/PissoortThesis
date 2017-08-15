#rsconnect::setAccountInfo(name, token, secret)  # FILL it !
## Connect to your rsconnect account
## Load the data if needed
# load(fit_stan.RData)

library(shinystan)
my_sso <- launch_shinystan(fit_stan) # Put you stan model in here.


## Put your shinyapps.io account, the name of your app.
deploy_shinystan(my_sso, appName = "ShinyStanGEV", account = "proto4426")




