library(PissoortThesis)

data("max_years")

library(ggplot2)
library(plotly)
library(gridExtra)
library(grid)
library(mgcv)


gam3.0 <- gam(Max ~ s(Year, k = 20), data = max_years$df, method = "REML")

# to generate random values from a multivariate normal
rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig) ;   m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

Vb <- vcov(gam3.0) # Bayesian covariance matrix of the model coefficients.
#it's conditional upon the smoothing parameter(s). add unconditional=T
#to adjust for the smoothing parameters being estimated rather than known values,


## Compute the coverages
'inCI' <- function(x, upr, lwr) {
  # =T if all evaluation points g lie within the interval and =F otherwise.
  all(x >= lwr & x <= upr)
}


shinyServer(function(input, output) {

  output$plot1 <- renderPlot({

    grid <- nrow(max_years$df) # Defines grid of values for which we will generate predictions.
    newd <- with(max_years$df,
                 data.frame(Year = seq(min(Year), max(Year), length = grid)))
    pred <- predict(gam3.0, newd, se.fit = TRUE)
    se.fit <- pred$se.fit


    set.seed(input$seed)
    # We want N draws from [\hat{beta}-beta, \hat{u}-u] which is ~multi N(0,Vb)
    BUdiff <- rmvn(input$sim, mu = rep(0, nrow(Vb)), sig = Vb)

    # Then compute \hat{f}(x)-f(x) which is  C_g%*%[\hat{beta}-beta, \hat{u}-u]
    Cg <- predict(gam3.0, newd, type = "lpmatrix")
    simDev <- Cg %*% t(BUdiff)

    # Find absolute values of the standardized deviations from the true model
    absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))

    # Max of the absolute standardized dev. at the grid of x values for each simul
    masd <- apply(absDev, 2L, max)

    # Find the crit value used to scale standard errors to yield the simultaneous interval
    crit <- quantile(masd, prob = 1-(input$level*.01), type = 8)

    # Now, compute and show the pointwise vs simultaneous confidence interval !
    pred <- transform(cbind(data.frame(pred), newd),
                      uprP = fit + (qnorm(1-input$level*.01/2) * se.fit),
                      lwrP = fit - (qnorm(1-input$level*.01/2) * se.fit),
                      uprS = fit + (crit * se.fit),
                      lwrS = fit - (crit * se.fit))


    sims <- rmvn(input$sim, mu = coef(gam3.0), sig = Vb)
    fits <- Cg %*% t(sims) # contains N draws from the posterior

    nrnd <- input$draws   ;    rnd <- sample(input$sim, nrnd)

    stackFits <- stack(as.data.frame(fits[, rnd]))
    stackFits <- transform(stackFits, Year = rep(newd$Year, length(rnd)))


    # Shows Coverages
    fitsInPCI <- apply(fits, 2L, inCI, upr = pred$uprP, lwr = pred$lwrP)
    fitsInSCI <- apply(fits, 2L, inCI, upr = pred$uprS, lwr = pred$lwrS)

    pointw_cov <- sum(fitsInPCI) / length(fitsInPCI)  # Point-wise
    simult_cov <- sum(fitsInSCI) / length(fitsInSCI)  # Simultaneous


    interval <- c("pointwise" =  "yellow", "simultaneous" = "darkred")

    ggplot(pred, aes(x = Year, y = fit)) +
      geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = "simultaneous"), alpha = 0.4) +
      geom_ribbon(aes(ymin = lwrP, ymax = uprP, fill = "pointwise"), alpha = 0.4) +
      geom_path(lwd = 2) +
      geom_path(data = stackFits, mapping = aes(y = values, x = Year, group = ind),
                alpha = 0.5, colour = "grey20") +
      labs(y = expression( Max~(T~degree*C)), x = "Year",
           title = "Point-wise & Simultaneous 95% conf. intervals for fitted GAM",
           subtitle = sprintf("Each line is one of %i draws from the posterior distribution of the model", input$draws)) +
      annotate(geom = "text", label = paste("coverages", " are : \n",
                                            round(pointw_cov, 5),
                                          " for pointwise \n", "   ",
                                          round(simult_cov, 5),
                                          " for simultaneous"),
               x = 1915,
               y = 33, col = "#33666C" , size = 6) +
      scale_fill_manual(name = "Interval", values = interval) +
      guides(colour = guide_legend(override.aes = list(size = 15))) +
      theme_piss(size_p = 25, plot.subtitle = text(31, hjust = 0.5, colour = "#33666C"),
                 legend.position = c(.888, .152))
    # expression(paste(underline("Coverage"), " is ", pointw_cov,
    #                  " for pointwise and ", simult_cov,
    #                  " for simultaneous")),
  })

})
#shinyApp(ui = ui, server = server)
