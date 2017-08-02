library(shiny)
library(shinydashboard)

#header <- dashboardHeader(div(style="color:#0000FF", p(title ="Extreme Value Test")))
header <- dashboardHeader( title = "Extreme Value")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("GEV distribution", tabName = "gev", icon = icon("dashboard"),
             badgeLabel = "chap.1", badgeColor = "green"),
    menuItem("Trend Models", icon = icon("thermometer-1"), tabName = "trend",
             badgeLabel = "chap.5", badgeColor = "green"),
    menuItem("Splines GAM", icon = icon("coffee"), tabName = "splines",
             badgeLabel = "chap.5", badgeColor = "red"),
    menuItem("Source code : Github", icon = icon("file-code-o"),
             href = "https://github.com/proto4426/PissoortThesis")
  )
)


body <- dashboardBody(
  tabItems(

   # First tab content
   tabItem(tabName = "gev",
           titlePanel(h2("Visualize the 3 GEV distributions")),

          fluidRow(column(3,
                          numericInput("ksi1", "Shape parameter for Weibull-type",
                                      "-0.5", min = "-100", max = "0", step = 1 ),
                          numericInput("ksi3", "Shape parameter for Fréchet-type",
                                      "0.5", min = "0", max = "100", step = 1 )
          ),
          width = "100px",

          column(3, offset = 1,
                 numericInput("mu", "Which location parameter ?",
                              "0", min = "-100000", max = "100000" ),
                 numericInput("sig", "Which scale parameter ?",
                              "1", min = "0", max = "10000" )
          ),

          mainPanel(
          plotOutput("plot1", height = '500px', width = "750px")
          )
         ), h2("Refer to Section 1.2.1 of the text for more information")),

  # Second tab content
  tabItem(tabName = "trend",
          titlePanel(h2("Visualize trend : first modelling")),
          fluidPage(selectInput("max",
                      "Do you want minima or maxima ? ", c('Max', 'Min') ),
          selectInput("fit",
                      label = "Which fitting method ? ",
                      choices = c("linear trend model" = "lm",
                                  "nonparam trend model" =  "loess",
                                  "broken linear trend" = "bl",
                                  "broken linear + linear trend" = "blll",
                                  "All 3 methods together" = 'all') ),

          mainPanel(
            plotOutput("plot2", height = '500px', width = "750px")
          )), h2("Refer to Section 5.2.2 of the text for more information")
          ),

  # Third tab content
  tabItem(tabName = "splines",
          titlePanel(h2("Simulate GAM fit for the trend to visualize uncertainty")),
          fluidRow(
            column(3,
                   numericInput("level", "Which level alpha? (in %) ",
                                "5", min = "0", max = "100" ),

                   numericInput("sim", "Howmuch Simulations M ? ",
                                "50", min = "2", max = "1000" )
            ),  width = "50px",
            fluidRow(
              column(3, offset = 1,
                     numericInput("seed", "Set the seed ",
                                  "99", min = "1", max = "1000000000" ),

                     numericInput("draws", "Draw simulations ? ( < M)",
                                  "50", min = "2", max = "1000"  )
              ),
              mainPanel(
                plotOutput("plot3", height = '500px', width = "750px")
              )
            )), h2("Refer to Section 5.2.3 of the text for more information")
      )
 )
)


ui <- dashboardPage(header, sidebar, body)


library(evd)
library(ggplot2)
library(PissoortThesis)

library(plotly)
library(gridExtra)
library(grid)

library(mgcv)

data("min_years")
data("max_years")


### For Third application (splines)

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



server <- function(input, output) {

  #### First application

  output$plot1 <- renderPlot({

    ## Create data frames for ggplot
    'GEVdfFun' <-
      function (x = seq(input$mu-10, input$mu + 10, length = 5e3), mu = 0, sig = 1, ksi = 0) {
        if (ksi ==0) dens <-  (sig^-1) * exp(-(x-mu)/sig) * exp(-exp(-(x-mu)/sig))
        else   s <- (1 + ksi * (x - mu)/sig)^(-(ksi)^-1 - 1)
        t <- (1 + ksi * (x - mu)/sig)^(-(ksi)^-1)
        if (ksi < 0) {dens <-  s * exp(-t) * ( (x - mu)/sig  < -1/ksi ) }
        if (ksi > 0) {dens <- sig^{-1} * s * exp(-t) * ( (x - mu)/sig  > -1/ksi ) }

        df <- data.frame(x = x, density = dens, xi = as.factor(ksi),
                         mu = as.factor(mu), scale = as.factor(sig))
        return(df)
      }


    ksi_gumb <- 0

    GEVdf <- rbind(GEVdfFun(mu = input$mu, sig=input$sig, ksi = input$ksi1),
                   GEVdfFun(mu = input$mu, sig=input$sig, ksi = ksi_gumb),
                   GEVdfFun(mu = input$mu, sig=input$sig, ksi = input$ksi3))

    # Dealing with endpoints
    endpoint_w <- input$mu - (input$sig / input$ksi1)
    endpoint_f <- input$mu - (input$sig / input$ksi3)

    dens_f <- ifelse(GEVdf[GEVdf$xi == input$ksi3,]$density < endpoint_f, NA,
                     GEVdf[GEVdf$xi == input$ksi3,]$density )
    GEVdf[GEVdf$xi == input$ksi3,]$density <- dens_f


    # plot the normal distribution as reference

    GEVdf <- cbind(GEVdf, norm = dnorm(GEVdf$x, mean = input$mu, sd = input$sig))

    #GEVdf[GEVdf$density < 10^{-312}, ]$density <- NA

    pres <- labs(title = expression(paste(underline(bold('Generalized Extreme Value density')))),
                 colour = expression(paste(xi,"=")))


    ggplot(GEVdf, aes(x = x, y = density, colour = xi )) +
      geom_line() + pres + theme_classic() +
      geom_line(aes(x = x, y = norm, col = "normal"), col = "black", linetype = 3) +
        theme(plot.title = element_text(size = 20, hjust=0.5,
                                        colour = "#33666C", face="bold"),
              axis.title = element_text(face = "bold", size= 15,
                                        colour = "#33666C"),
              axis.line = element_line(color="#33666C", size = .45),
              legend.background = element_rect(colour = "black"),
              legend.key = element_rect(fill = "white") ) +
      theme(legend.title = element_text(colour="#33666C",
                                        size=18, face="bold")) +
      theme(legend.key = element_rect(colour = "black")) +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      geom_point(aes(x = endpoint_f, y = 0),size = 3.5) +
      geom_point(aes(x = endpoint_w, y = 0), col="red",size = 3.5)

  })


  #### Second application

  output$plot2 <- renderPlot({
    x <- cbind.data.frame(max_years$df, Min = min_years$data)


    g <- ggplot(x, aes_string(x = 'Year', y = input$max)) +
      theme_bw() + geom_line() + geom_point() +
      theme(plot.title = element_text(size = 22, hjust=0.5,
                                      colour = "#33666C", face="bold"),
            axis.title = element_text(face = "bold", size= 18,
                                      colour = "#33666C"),
            axis.line = element_line(color="#33666C", size = .45),
            legend.position = c(.888, .152),
            legend.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = "white") ) +
      scale_colour_manual(name="Trend",
                          values=c(Linear="blue", BrokenLinear="cyan", LOESS="red")) +
      guides(colour = guide_legend(override.aes = list(size = 4)))


    g_linear <- geom_smooth(method='lm',formula=y~x, aes(colour = "Linear"))

    g_smooth <-  stat_smooth(method = "loess", se = F, aes(colour = 'LOESS'))

    g_bl_max <-  list(geom_line(data = max_years$df[max_years$df$Year %in% 1901:1975,],
                                aes(x = Year, colour = "BrokenLinear",
                                    y = predict(lm(max_years$data[1:75] ~ max_years$df$Year[1:75]))),
                                size = 1.5, linetype = "twodash"),
                      geom_line(data = max_years$df[max_years$df$Year %in% 1977:2016,],
                                aes(x = Year, colour = "BrokenLinear",
                                    y = predict(lm(max_years$data[77:116] ~ max_years$df$Year[77:116]))),
                                size = 1.5, linetype = "twodash")     )

    g_bl_min <- list(geom_line(data = min_years$df[min_years$df$Year %in% 1901:1975,],
                               aes(x = Year, colour = "BrokenLinear",
                                   y = predict(lm(min_years$data[1:75] ~ min_years$df$Year[1:75]))),
                               size = 1.5, linetype = "twodash"),
                     geom_line(data = min_years$df[min_years$df$Year %in% 1977:2016,],
                               aes(x = Year, colour = "BrokenLinear",
                                   y = predict(lm(min_years$data[77:116] ~ min_years$df$Year[77:116]))),
                               size = 1.5, linetype = "twodash")     )


    if(input$fit == 'lm')  g <- g + g_linear
    if(input$fit == 'loess') g <- g + g_smooth


    ## Broken linear trends

    if (input$max == "Max" & input$fit == "bl")
      g <- g + g_bl_max

    if (input$max == "Min" & input$fit == "bl")
      g <- g + g_bl_min


    ## Broken linear and linear trend

    if (input$max == "Max" & input$fit == "blll")
      g <- g + g_linear + g_bl_max

    if (input$max == "Min" & input$fit == "blll")
      g <- g + g_linear + g_bl_min


    ### ALL the methods

    if (input$max == "Max" & input$fit == "all")
      g <- g + g_linear + g_bl_max + g_smooth

    if (input$max == "Min" & input$fit == "all")
      g <- g + g_linear +  g_bl_min + g_smooth



    if(input$max == "Min") g <- g +
      labs(title = "Complete Serie of Annual TN in Uccle")
    else g <- g + labs(title = "Complete Serie of Annual TX in Uccle")
    g
  })


  #### Third application

  output$plot3 <- renderPlot({

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
           subtitle = sprintf("Each line is one of %i draws to display (randomly) from the posterior distribution of the model", input$draws)) +
      annotate(geom = "text", label = paste("coverages", " are : \n",
                                            round(pointw_cov, 5),
                                            " for pointwise \n", "   ",
                                            round(simult_cov, 5),
                                            " for simultaneous"),
               x = 1915,
               y = 33, col = "#33666C" , size = 6) +
      scale_fill_manual(name = "Interval", values = interval) +
      guides(colour = guide_legend(override.aes = list(size = 15))) +
      theme_light() +
      theme(plot.title = element_text(size = 22, hjust=0.5,
                                      colour = "#33666C", face="bold"),
            axis.title = element_text(face = "bold", size= 18,
                                      colour = "#33666C"),
            axis.line = element_line(color="#33666C", size = .45),
            legend.position = c(.888, .152),
            legend.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = "white"),
            plot.subtitle = text(31, hjust = 0.5,
                                 colour = "#33666C")
            )

    # expression(paste(underline("Coverage"), " is ", pointw_cov,
    #                  " for pointwise and ", simult_cov,
    #                  " for simultaneous")),
  })

}

shinyApp(ui, server)



