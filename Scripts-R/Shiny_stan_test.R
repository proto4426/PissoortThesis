#rsconnect::setAccountInfo("proto4426", , secret)  # FILL it !
## Connect to your rsconnect account
## Load the data if needed
# load(fit_stan.RData)
library(rsconnect)
setAccountInfo(name='proto4426', token='60D7FBC0924BAEDBD1B19F21294F8CF2', secret='yJyLnAmErkzIYLOlHlIuWa+bQrvVJovkEP8nXBuj')

library(shinystan)
gev_app <- launch_shinystan(fit_stan) # Put you stan model in here.


## Put your shinyapps.io account, the name of your app.
deploy_shinystan(gev_app, appName = "ShinyStanGEV_converged",
                 account = "proto4426",deploy = T)



library(RJSONIO)
library(pkgdll)
library(rjsonioUser)
