library(htmltools)
library(shinythemes)
library(shiny)
library(shinydashboard)


header <- dashboardHeader( title = h4("Extreme Values' Analysis in Uccle :", strong("PisssortThesis")),
                           titleWidth = 450)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("GEV distribution", tabName = "gev", icon = icon("dashboard"),
             badgeLabel = "chap.1", badgeColor = "green"),
    menuItem("Trend Models", icon = icon("thermometer-1"), tabName = "trend",
             badgeLabel = "chap.5", badgeColor = "green"),
    menuItem("Splines GAM", icon = icon("coffee"), tabName = "splines",
             badgeLabel = "chap.5", badgeColor = "red"),
    menuItem("Neural Networks", icon = icon("plane"), tabName = "NN",
             badgeLabel = "chap.6", badgeColor = "red"),
    menuItem("Bayesian", icon = icon("plane"), tabName = "bay",
             badgeLabel = "chap.7", badgeColor = "blue"),
    menuItem("Source code : Github", icon = icon("file-code-o"),
             href = "https://github.com/proto4426/PissoortThesis")
  )
)

# 'tabItems' <- function (...)  {
#   lapply(..., tagAssert)
#   div(class = "tab-content", ...)
# }

body <- dashboardBody(
 includeCSS("custom/style.css"),
  tabItems(

   # First tab content
   tabItem(tabName = "gev",
           titlePanel(h2("Visualize the 3 GEV distributions")),

          fluidRow(column(3,
                          numericInput("ksi1", "Shape parameter for Weibull-type",
                                      "-0.5", min = "-100", max = "0", step = .1 ),
                          numericInput("ksi3", "Shape parameter for Fréchet-type",
                                      "0.5", min = "0", max = "100", step = .1 )
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
      ),

  # Fourth tab content
  tabItem(tabName = "NN",
   # navbarPage(#theme = shinytheme("united"),
   # "EVT thesis of Antoine Pissoort",
     #tabPanel("GEV-CDN : Neural Networks to fit flexible nonstationary Nonlinear GEV Models in parallel + prevent overfitting (bagging)",
   titlePanel(h3("GEV-CDN : Neural Networks to fit flexible nonstationary Nonlinear GEV Models in parallel + prevent overfitting (bagging)")),
      sidebarPanel(
       fluidRow(
        wellPanel(h3("Hyperparameters "),
                  selectInput("param",
                              label = "Which model ? ",
                              choices = c("Stationary" = "sta",
                                          "nonstationary in Location" =  "lin1",
                                          "nonstationary in Location and Scale" = "lin2",
                                          "nonstationary in Location, Scale and Shape" = "lin3") ),
                  selectInput("m",
                              label = "Activation function ? ",
                              choices = c("Identity" = "identity",
                                          "Logistic sigmoid" = "logsig",
                                          "Hyperbolic tangent" =  "tanh")),
                  code("Problem could occur for hyperbolic tan."),
                  sliderInput("hidd", "Number of hidden layers ?",
                              value = 0, min =0, max = 5 )
        ),
        wellPanel(numericInput("regprior",
                               "Std deviation for Gaussian prior on the weights ? (regularization)" ,
                               value = 1e5, min = 1e-5, max = 1e15, step = 10  ),
                  numericInput("beta1", "Parameters for the shifted Beta prior distribution (1) ?" ,
                               "6", min = "2", max = "1000"  ),
                  numericInput("beta2", "Parameters for the shifted Beta prior distribution (2) ?" ,
                               "9", min = "2", max = "1000"  )
        ),
        htmltools::p(actionButton("runfit","RUN GEV-CDN", icon("random") ),
                     align = "center", width = 9 ),
        sliderInput("bag", "Bagging resamples B ? ",
                    value = 1, min = 1, max = 500, step = 4 ),
        h5("1=no bagging. Only for nonstationary models only with hidden layer >0"),
        wellPanel(h3("Bootstrap conf. intervals "),
                  sliderInput("nboot", "Bootstrap resamples?",
                              value = 1, min = 1, max = 500, step = 4 ),
                  selectInput("method",
                              label = "Bootstrap method ?",
                              choices = c("Residual" = "residual",
                                          "Parametric" =  "parametric")),
                  htmltools::p(actionButton("runboot","RUN Bootstrap intervals", icon("random") ),
                               align = "center", width = 9 ),

                  checkboxInput("comp", "Comparison of the resdiual and parametric methods ? ", FALSE)
        )
      )
     ),
   mainPanel(
     fluidRow(
     tabsetPanel(#h3( "If bootstrap, look at console to follow computation's progress "),
       tabPanel(h5("Model fitting + parameters' intervals"),
                plotOutput("plotfitNN", height = '500px', width = "800px"),
                br(),        br(),
                h4(strong("Computation times (sec.)")),
                tableOutput("datatable"),

                plotOutput("plotNNboot", height = '700px', width = "600px")
       ),
       tabPanel("Interval's comparison + Summary Table",
                h4(strong("Difference between residual and parametric Bootrstrap quantiles")),
     plotOutput("plotBootcomp", height = '300px', width = "800px"),
     br(),        br(),
     h5(strong("Parameters' Estimated values from the selected model")),
     DT::dataTableOutput("datatableFin", width = "800px")
   ),
   tabPanel("Informations", icon = icon("info-circle"),
            htmlOutput("infoNN")  # See start of server()
   )
     )
     )
   )
  ),
  tabItem(tabName = "bay",
          titlePanel(h3("Bayesian Analysis")),
                   sidebarPanel(
                     wellPanel(h2("Model"),
                               sliderInput("iterchain", "Number of iterations by chains ",
                                           value = 500, min = 10, max = 1e4, step = 20),
                               numericInput("seed", "Set the seed ",
                                            value = 123, min = 1, max = 1e15)
                     ),
                     wellPanel(h2("Priors"),
                               fluidRow(column(3,
                                               numericInput("priormumean", "mean mu",
                                                            value = 30, min = -100, max = 200, step = 1),
                                               numericInput("priormusd", "SD mu",
                                                            value = 40, min = 1e-15, max = 1e15, step = 2)
                               ), column(3,
                                         numericInput("priormu1mean", "mean mu1",
                                                      value = 0, min = -1000, max = 1000, step = 1),
                                         numericInput("priormu1sd", "SD for mu1",
                                                      value = 40, min = 1e-15, max = 5e3, step = 2)
                               ), column(3,
                                         # ), column(6,
                                         numericInput("priorlogsigmean", "mean logsig",
                                                      value = 0, min = 1, max = 10, step = 1),
                                         numericInput("priorlogsigsd", "SD logsig",
                                                      value = 10, min = 1e-15, max = 5e3, step = 2)
                               ),
                               column(3,
                                      numericInput("priorximean", "mean xi",
                                                   value = 0, min = -100, max = 100, step = 0.1),
                                      numericInput("priorxisd", "SD xi",
                                                   value = 10, min = 1e-15, max = 5e3, step = 2)
                               )),
                               h5("Take them uninformative by default ")
                     ),
                     htmltools::p(actionButton("run","RUN Gibbs Sampler", icon("random") ),
                                  align = "center", width = 9 ),
                     wellPanel(h2("Diagnostics"),
                               sliderInput("start", "Number of chains with different starting values ?",
                                           value = 2, min = 1, max = 10, step = 1),
                               sliderInput("burnin", "Number of burnin by chains ",
                                           value = 10, min = 0, max = 5e3, step = 10),
                               column(4,
                                      checkboxInput("gelman", "Gelman-R", F),
                                      checkboxInput("geweke", "Geweke", F)
                               ),
                               column(4,
                                      checkboxInput("autocor", "Autcorr", F),
                                      checkboxInput("crosscor", "Cross-corr", F)
                               ),
                               column(4, checkboxInput("raft", "Raftery-Coda", F)
                               )
                     ),
                     br(),
                     wellPanel(h2("Posterior Predictive"),
                               sliderInput("from", "Start at year ",
                                           value = 1901, min = 1901, max = 2016, step = 1 ),
                               sliderInput("fut", "Prediction in the future ? ",
                                           value = 1, min = 0, max = 250, step = 5 ),
                               numericInput("dens", "Show densities every which years ?",
                                            value = "10", min = "1", max = "30" ),
                               checkboxInput("show", "Show the intervals' lengths on the graph ? ", FALSE)
                     )
                   ),
                   mainPanel(
                     tabsetPanel(
                       tabPanel(h4("Predictive Posterior"),
                                br(),
                                wellPanel(h5("This application demonstrates the Bayesian Results. Information are given on the ",
                                             icon("info-circle"), "Informations tab."),
                                          h5(strong("Click on the"), icon("random"), strong("RUN button to compute the results"))
                                ),
                                br(),        br(),
                                plotOutput("plotPred1", height = '800px', width = "800px"),
                                br(),        br(),
                                plotOutput("plotPred2", height = '500px', width = "800px"),

                                verbatimTextOutput("text")
                       ), # tabPanel
                       tabPanel(h4("MCMC Diagnostics"),
                                br(),
                                plotOutput("plot.chains", height = '300px', width = "800px"),
                                plotOutput("plot.chains2", height = '300px', width = "800px"),
                                h4(strong("Starting values : ")),
                                DT::dataTableOutput("DTstart", width = "750px"),
                                h4(strong("Acceptance rates : ")),
                                DT::dataTableOutput("accrates", width = "600px"),
                                code("These values are recommended to be around 0.4 (chosen automatically here). If this is not the case, convergence is expected to be slower."),
                                br(), br(),   # ==================================================================
                                htmltools::p(h4(strong("Other MCMC Diagnostics")), align = "center"),
                                tabsetPanel(#"Other MCMC Diagnostics",
                                  tabPanel("Gelamn-Rubin",
                                           plotOutput("gelman", height = "500px", width = "750px")
                                  ),
                                  tabPanel("Correlation",
                                           plotOutput(outputId="autocorr", width="600px",height="400px"),
                                           plotOutput(outputId="crosscorr", width="450px",height="400px")
                                  ),
                                  tabPanel("Geweke",
                                           plotOutput(outputId="geweke", width="650px",height="500px")
                                  ),
                                  tabPanel("Raftery-Coda",
                                           column(6,
                                                  DT::dataTableOutput("raft", width = "400px")),
                                           column(6, htmlOutput("info_raft")
                                           )
                                  )
                                )
                       ),
                       tabPanel("Informations", icon = icon("info-circle"),
                                htmlOutput("infobay")#,
                       )
                     )  # tabsetPanel
                   )  # mainPanel
          ) # tabPanel
  )
)
 # )
#)



ui <- dashboardPage(header, sidebar, body)


library(evd)
#devtools::install_github("proto4426/PissoortThesis", force = T)
library(PissoortThesis)
library(plotly)
library(gridExtra)
library(grid)

library(GEVcdn)
library(foreach)
library(doParallel)
library(doSNOW)
library(tidyverse)
library(gridExtra)
library(grid)
library(DT)
library(pander)
library(mgcv)

library(reshape2)
library(coda)
library(mvtnorm)
library(HDInterval)
library(ggjoy)
library(viridis)
library(ggcorrplot)


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



# NN : tanh activation function
'hyp.tan' <- function(x) ( exp(x) - exp(x) ) / ( exp(x) + exp(-x) )

# Create a function to easily compare all the methods. df contains the residuals for all quantiles
# See output$plot3
'compare_methods_quantiles_plot' <- function(df) {
  col.quantiles <- c("97.5%" = "red", "50%" = "black", "2.5%" = "blue")
  gloc.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
    geom_line(aes(y = loc.2.5., col  = "2.5%")) +
    geom_line(aes(y = loc.97.5., col  = "97.5%"))+
    labs(title = "Location parameter") +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    theme_piss()
  gsc.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
    geom_line(aes(y = scale.2.5., col  = "2.5%")) +
    geom_line(aes(y = scale.97.5., col  = "97.5%")) +
    labs(title = "Scale parameter") +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    theme_piss()
  gsh.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
    geom_line(aes(y = shape.2.5., col  = "2.5%")) +
    geom_line(aes(y = shape.97.5., col  = "97.5%")) +
    labs(title = "Shape parameter") +
    scale_colour_manual(name = "Quantiles", values = col.quantiles) +
    theme_piss()
  PissoortThesis::grid_arrange_legend(gloc.lin, gsc.lin, gsh.lin)
}


server <- function(input, output) {
  #browser()

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


  #### Third application (Splines)

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
  })


  #### Fourth application (NN)

  fit <- eventReactive(input$runfit, {

    if(input$m %in% c("tanh", "logsig") )
      validate(
        need(input$hidd > 0,
             label = "Activation function is nonlinear --> hidden layers must be >0.")
      )

    if(input$m == "identity")  Th <- gevcdn.identity
    else if(input$m == "logsig")  Th <- gevcdn.logistic
    else if(input$m == "tanh")  Th <- hyp.tan

    if(input$param == "sta")  fixed <- c("location", "scale", "shape")
    else if(input$param == "lin1")  fixed <- c("scale", "shape")
    else if(input$param == "lin2")  fixed <- c("shape")
    else if(input$param == "lin3")  fixed <- NULL


    x <- as.matrix(1:(length(max_years$data)))
    y <- as.matrix(max_years$data)

    if(input$bag == 1 ) {

      t <- proc.time()
      weights <- PissoortThesis::gevcdn.fit2(x = x, y = y,
                                             n.trials = 1,
                                             n.hidden = input$hidd,
                                             Th = Th ,
                                             fixed = fixed ,
                                             beta.p = input$beta1,
                                             beta.q = input$beta2,
                                             sd.norm = input$regprior,
                                             silent = T)
      time.fit <- (proc.time()-t)[3]



      ### Compute the quantiles an plot the results
      parms.best <- gevcdn.evaluate(x, weights)

      q.best <- sapply(c(0.025, 0.05, 0.1,  0.5, 0.9, 0.95, 0.975), VGAM::qgev,
                       location = parms.best[,"location"],
                       scale = parms.best[,"scale"],
                       shape = parms.best[,"shape"])

      df <- data.frame(year = x , obs = y, q.025 = q.best[,1], q.05 = q.best[,2],
                       q.10 = q.best[,3],
                       q.50 = q.best[,4], q.90 = q.best[,5], q.975 = q.best[,7])


    }  else {

      withProgress(message = 'Baggig computation', value = 0, {

        'mean_of_list' <- function(param = parms.on) {
          parms <- matrix(0, nrow = nrow(param[[1]]), ncol = 3)
          for (i in 1:length(param)){
            parms <- parms + as.matrix(param[[i]])
          }
          parms <- parms / length(param)
          parms
        }

        M <- input$bag
        H <- input$hidd
        sdnorm <- input$regprior


        cores <- detectCores()
        cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
        registerDoSNOW(cl)
        pb <- txtProgressBar(max = M, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)

        t <- proc.time()

        incProgress(0.5)

        bag_par <- foreach(i = 1:M,
                           .packages = c("PissoortThesis", "GEVcdn"),
                           .options.snow = opts) %dopar% {

                             weights.on <-  PissoortThesis::gevcdn.bag2(x = x, y = y,
                                                                        iter.max = 100,
                                                                        fixed = fixed,
                                                                        #iter.step = 10,
                                                                        n.bootstrap = 2,
                                                                        n.hidden = H,
                                                                        sd.norm = sdnorm,
                                                                        Th = Th,
                                                                        silent = F)
                             parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)

                             mean_of_list(parms.on)
                           }
        close(pb)
        stopCluster(cl)
        time.fit <- unname(as.vector((proc.time()-t)))[3]

        setProgress(1)


      })

      bag_par_on <- mean_of_list(bag_par)


      q <- t(sapply(y, quantile, probs = c(.1, .5, .9)))

      q.025.on <- q.05.on <- q.10.on <- q.50.on <- q.90.on <- q.95.on <- q.975.on <- c()
      for(i in seq_along(bag_par)){
        q.025.on <- cbind(q.025.on, VGAM::qgev(p = 0.025,
                                               location = bag_par[[i]][,"location"],
                                               scale = bag_par[[i]][,"scale"],
                                               shape = bag_par[[i]][,"shape"]))
        q.05.on <- cbind(q.05.on, VGAM::qgev(p = 0.05,
                                             location = bag_par[[i]][,"location"],
                                             scale = bag_par[[i]][,"scale"],
                                             shape = bag_par[[i]][,"shape"]))
        q.10.on <- cbind(q.10.on, VGAM::qgev(p = 0.1,
                                             location = bag_par[[i]][,"location"],
                                             scale = bag_par[[i]][,"scale"],
                                             shape = bag_par[[i]][,"shape"]))
        q.50.on <- cbind(q.50.on, VGAM::qgev(p = 0.5,
                                             location = bag_par[[i]][,"location"],
                                             scale = bag_par[[i]][,"scale"],
                                             shape = bag_par[[i]][,"shape"]))
        q.90.on <- cbind(q.90.on, VGAM::qgev(p = 0.9,
                                             location = bag_par[[i]][,"location"],
                                             scale = bag_par[[i]][,"scale"],
                                             shape = bag_par[[i]][,"shape"]))
        q.95.on <- cbind(q.975.on, VGAM::qgev(p = 0.95,
                                              location = bag_par[[i]][,"location"],
                                              scale = bag_par[[i]][,"scale"],
                                              shape = bag_par[[i]][,"shape"]))
        q.975.on <- cbind(q.975.on, VGAM::qgev(p = 0.975,
                                               location = bag_par[[i]][,"location"],
                                               scale = bag_par[[i]][,"scale"],
                                               shape = bag_par[[i]][,"shape"]))
      }

      df <- data.frame(year = 1901:2016, obs = y, q.025 = rowMeans(q.025.on),
                       q.10 = rowMeans(q.10.on), q.50 = rowMeans(q.50.on),
                       q.90 = rowMeans(q.90.on), q.975 = rowMeans(q.975.on))
    }

    list(df = df, time.fit = time.fit)

  })


  output$plotfitNN <- renderPlot({

    df <- fit()[["df"]]

    # in ggplot

    col.quantiles <- c("2.5% and 97.5%" = "blue", "10% and 90%" = "green", "50%" = "red")
    gg.cdn <- ggplot(df, aes(x = year, y = obs)) +
      geom_line() + geom_point() +
      geom_line(aes(y = q.025, col = "2.5% and 97.5%")) +
      geom_line(aes(y = q.50, col = "50%")) +
      geom_line(aes(y = q.975, col = "2.5% and 97.5%")) +
      geom_line(aes(y = q.10, col = "10% and 90%")) +
      geom_line(aes(y = q.90, col = "10% and 90%")) +
      scale_colour_manual(name = "Quantiles", values = col.quantiles) +
      labs(title = "GEV-CDN quantiles with identity link for the location µ(t)",
           y =  expression( Max~(T~degree*C))) +
      theme_piss()
    gg.cdn

  })

  boot <- eventReactive( input$runboot, {

    withProgress(message = 'Bootstrap computation', value = 0, {

      if(input$m == "identity")  Th <- gevcdn.identity
      else if(input$m == "logsig")  Th <- gevcdn.logistic
      else if(input$m == "tanh")  Th <- hyp.tan

      if(input$param == "sta")   fixed <- c("location", "scale", "shape")
      else if(input$param == "lin1")   fixed <- c("scale", "shape")
      else if(input$param == "lin2")   fixed <- c("shape")
      else if(input$param == "lin3")   fixed <- NULL


      x <- as.matrix(1:(length(max_years$data)))
      y <- as.matrix(max_years$data)

      if(input$method == "residual") meth = "residual"
      else if (input$method == "parametric")  meth <- "parametric"


      ###### Bootstrap confidence intervals
      B <- input$nboot

      H <- input$hidd

      cores <- detectCores()
      cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = B, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)

      t <- proc.time()

      incProgress(0.25)

      boot_par <- foreach::foreach(i = 1:B,
                                   .packages = c("PissoortThesis", "GEVcdn"),
                                   .options.snow = opts) %dopar% {
                                     ci.lin <- gevcdn.bootstrap(n.bootstrap = 2,
                                                                x = x, y = y,
                                                                iter.max = 100,
                                                                n.hidden = H,
                                                                fixed = fixed,
                                                                Th = Th,
                                                                n.trials = 1,
                                                                boot.method = meth)
                                     list(loc = ci.lin$location.bootstrap,
                                          sc = ci.lin$scale.bootstrap,
                                          sh = ci.lin$shape.bootstrap
                                     )
                                   }
      close(pb)
      parallel::stopCluster(cl)
      time.boot <- unname(as.vector((proc.time()-t)))[3]
      incProgress(0.75)

      ## Aggregate the results of the resamples and the parralel computation

      boot.loc.h <- matrix(nrow = nrow(boot_par[[1]]$loc))
      for (i in 1:B)   boot.loc.h <- cbind(boot.loc.h, boot_par[[i]]$loc)
      boot.loc.h_f <- boot.loc.h[,-1]

      boot.sc.h <- matrix(nrow = nrow(boot_par[[1]]$loc))
      for (i in 1:B)   boot.sc.h <- cbind(boot.sc.h, boot_par[[i]]$sc)
      boot.sc.h_f <- boot.sc.h[,-1]

      boot.sh.h <- matrix(nrow = nrow(boot_par[[1]]$loc))
      for (i in 1:B)   boot.sh.h <- cbind(boot.sh.h, boot_par[[i]]$sh)
      boot.sh.h_f <- boot.sh.h[,-1]


      ## Compute the quantiles
      b.loc.h <- t(apply(boot.loc.h_f, 1, quantile, p = c(0.025, 0.5, 0.975)))
      b.sc.h <- t(apply(boot.sc.h_f, 1, quantile, p = c(0.025, 0.5, 0.975)))
      b.sh.h <- t(apply(boot.sh.h_f, 1, quantile, p = c(0.025, 0.5, 0.975)))

      setProgress(1)

    })

    df.boot <- data.frame(year = 1901:2016 , obs = y,
                          loc = b.loc.h,
                          scale = b.sc.h,
                          shape = b.sh.h)

    list(df.boot = df.boot, time.boot = time.boot, meth = meth)

  })



  output$plotNNboot <- renderPlot({

    meth <- boot()[["meth"]]
    df.boott <- boot()[["df.boot"]]


    # for the location
    gg.boot.loc <- ggplot(df.boott, aes(x = year)) +
      geom_line(aes(y = loc.2.5.), col = "blue") +
      geom_line(aes(y = loc.97.5.), col = "blue") +
      labs(title = "For the Location parameter", y = "") +
      theme_piss()
    # for the scale
    gg.boot.sc <- ggplot(df.boott, aes(x = year)) +
      geom_line(aes(y = scale.2.5.), col = "green") +
      geom_line(aes(y = scale.97.5.), col = "green") +
      labs(title = "For the Scale parameter", y = "") +
      theme_piss()
    # for the shape
    gg.boot.sh <- ggplot(df.boott, aes(x = year)) +
      geom_line(aes(y = shape.2.5.), col = "red") +
      geom_line(aes(y = shape.97.5.), col = "red") +
      labs(title = "For the Shape parameter", y = "") +
      theme_piss()

    # All in one for the identity link model on the location
    grid.arrange(gg.boot.loc, gg.boot.sc, gg.boot.sh, ncol = 1,
                 top = textGrob(sprintf("GEV-CDN intervals with %s bootstrap", meth),
                                gp = gpar(fontsize = 24, font = 4, col ="black")))


  })


  output$plotBootcomp <- renderPlot({

    validate(
      need(input$comp,  "Check the 'comparison...' box to see the graphs")
    )

    if (input$comp) {

      if(input$method == "residual") {

        withProgress(message = 'Bootstrap 2nd computation', value = 0, {
          incProgress(0.25)


          ## Quantitative comparisons between the residual and parametric bootstraps

          df.boot.res <- boot()[["df.boot"]]


          ## Compute the parametric bootstrap

          # ===================================================================================

          if(input$m == "identity")  Th <- gevcdn.identity
          else if(input$m == "logsig")  Th <- gevcdn.logistic
          else if(input$m == "tanh")  Th <- hyp.tan

          if(input$param == "sta")   fixed <- c("location", "scale", "shape")
          else if(input$param == "lin1")   fixed <- c("scale", "shape")
          else if(input$param == "lin2")   fixed <- c("shape")
          else if(input$param == "lin3")   fixed <- NULL

          x <- as.matrix(1:(length(max_years$data)))
          y <- as.matrix(max_years$data)

          B <- input$nboot
          H <- input$hidd
          cores <- detectCores()
          # cl <- makeCluster(input$cores)
          cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
          #registerDoParallel(cl)
          registerDoSNOW(cl)
          pb <- txtProgressBar(max = B, style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)

          incProgress(0.5)

          boot_par_comp <- foreach::foreach(i = 1:B,
                                            .packages = c("PissoortThesis", "GEVcdn"),
                                            .options.snow = opts) %dopar% {
                                              ci.lin <- gevcdn.bootstrap(n.bootstrap = 2,
                                                                         x = x, y = y,
                                                                         iter.max = 100,
                                                                         n.hidden = H,
                                                                         fixed = fixed,
                                                                         Th = Th,
                                                                         n.trials = 1,
                                                                         boot.method = "parametric")
                                              list(loc = ci.lin$location.bootstrap,
                                                   sc = ci.lin$scale.bootstrap,
                                                   sh = ci.lin$shape.bootstrap
                                              )
                                            }
          close(pb)
          parallel::stopCluster(cl)

          incProgress(0.75)

          ## Aggregate the results of the resamples and the parralel computation

          boot.loc.par <- matrix(nrow = nrow(boot_par_comp[[1]]$loc))
          for (i in 1:B)   boot.loc.par <- cbind(boot.loc.par, boot_par_comp[[i]]$loc)
          boot.loc.par_f <- boot.loc.par[,-1]

          boot.sc.par <- matrix(nrow = nrow(boot_par_comp[[1]]$loc))
          for (i in 1:B)   boot.sc.par <- cbind(boot.sc.par, boot_par_comp[[i]]$sc)
          boot.sc.par_f <- boot.sc.par[,-1]

          boot.sh.par <- matrix(nrow = nrow(boot_par_comp[[1]]$loc))
          for (i in 1:B)   boot.sh.par <- cbind(boot.sh.par, boot_par_comp[[i]]$sh)
          boot.sh.par_f <- boot.sh.par[,-1]


          ## Compute the quantiles
          b.loc.h <- t(apply(boot.loc.par_f, 1, quantile, p = c(0.025, 0.5, 0.975)))
          b.sc.h <- t(apply(boot.sc.par_f, 1, quantile, p = c(0.025, 0.5, 0.975)))
          b.sh.h <- t(apply(boot.sh.par_f, 1, quantile, p = c(0.025, 0.5, 0.975)))


          df.boot.par <- data.frame(year = 1901:2016 , obs = y,
                                    loc = b.loc.h,
                                    scale = b.sc.h,
                                    shape = b.sh.h)
          # ===================================================================================


          diff.boot.compp <- df.boot.res[,-(1:2)] - df.boot.par[, -(1:2)]

          setProgress(1)

        })

        compare_methods_quantiles_plot(diff.boot.compp)

      }

      else if(input$method == "parametric"){

        withProgress(message = 'Bootstrap 2nd computation', value = 0, {
          incProgress(0.25)


          ## Quantitative comparisons between the residual and parametric bootstraps

          df.boot.res <- boot()[["df.boot"]]

          ## Compute the parametric bootstrap

          # ===================================================================================

          if(input$m == "identity")  Th <- gevcdn.identity
          else if(input$m == "logsig")  Th <- gevcdn.logistic
          else if(input$m == "tanh")  Th <- hyp.tan

          if(input$param == "sta")   fixed <- c("location", "scale", "shape")
          else if(input$param == "lin1")   fixed <- c("scale", "shape")
          else if(input$param == "lin2")   fixed <- c("shape")
          else if(input$param == "lin3")   fixed <- NULL

          x <- as.matrix(1:(length(max_years$data)))
          y <- as.matrix(max_years$data)

          B <- input$nboot
          H <- input$hidd
          cores <- detectCores()
          cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
          registerDoSNOW(cl)
          pb <- txtProgressBar(max = B, style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)

          incProgress(0.5)

          boot_par_comp <- foreach::foreach(i = 1:B,
                                            .packages = c("PissoortThesis", "GEVcdn"),
                                            .options.snow = opts) %dopar% {
                                              ci.lin <- gevcdn.bootstrap(n.bootstrap = 2,
                                                                         x = x, y = y,
                                                                         iter.max = 100,
                                                                         n.hidden = H,
                                                                         fixed = fixed,
                                                                         Th = Th,
                                                                         n.trials = 1,
                                                                         boot.method = "parametric")
                                              list(loc = ci.lin$location.bootstrap,
                                                   sc = ci.lin$scale.bootstrap,
                                                   sh = ci.lin$shape.bootstrap
                                              )
                                            }
          close(pb)
          parallel::stopCluster(cl)

          ## Aggregate the results of the resamples and the parralel computation

          boot.loc.par <- matrix(nrow = nrow(boot_par_comp[[1]]$loc))
          for (i in 1:B)   boot.loc.par <- cbind(boot.loc.par, boot_par_comp[[i]]$loc)
          boot.loc.par_f <- boot.loc.par[,-1]

          boot.sc.par <- matrix(nrow = nrow(boot_par_comp[[1]]$loc))
          for (i in 1:B)   boot.sc.par <- cbind(boot.sc.par, boot_par_comp[[i]]$sc)
          boot.sc.par_f <- boot.sc.par[,-1]

          boot.sh.par <- matrix(nrow = nrow(boot_par_comp[[1]]$loc))
          for (i in 1:B)   boot.sh.par <- cbind(boot.sh.par, boot_par_comp[[i]]$sh)
          boot.sh.par_f <- boot.sh.par[,-1]


          ## Compute the quantiles
          b.loc.h <- t(apply(boot.loc.par_f, 1, quantile, p = c(0.025, 0.5, 0.975)))
          b.sc.h <- t(apply(boot.sc.par_f, 1, quantile, p = c(0.025, 0.5, 0.975)))
          b.sh.h <- t(apply(boot.sh.par_f, 1, quantile, p = c(0.025, 0.5, 0.975)))


          df.boot.par <- data.frame(year = 1901:2016 , obs = y,
                                    loc = b.loc.h,
                                    scale = b.sc.h,
                                    shape = b.sh.h)
          # ===================================================================================


          diff.boot.compp <- df.boot.res[,-(1:2)] - df.boot.par[, -(1:2)]


          incProgress(0.75)

          # Create a function to easily compare all the methods. df contains the residuals for all quantiles
          'compare_methods_quantiles_plot' <- function(df) {
            col.quantiles <- c("97.5%" = "red", "50%" = "black", "2.5%" = "blue")
            gloc.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
              geom_line(aes(y = loc.2.5., col  = "2.5%")) +
              geom_line(aes(y = loc.97.5., col  = "97.5%")) +
              labs(title = "Location parameter") +
              scale_colour_manual(name = "Quantiles", values = col.quantiles) +
              theme_piss()
            gsc.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
              geom_line(aes(y = scale.2.5., col  = "2.5%")) +
              geom_line(aes(y = scale.97.5., col  = "97.5%")) +
              labs(title = "Scale parameter") +
              scale_colour_manual(name = "Quantiles", values = col.quantiles) +
              theme_piss()
            gsh.lin <- ggplot(cbind.data.frame(df, Year = 1901:2016), aes(x = Year)) +
              geom_line(aes(y = shape.2.5., col  = "2.5%")) +
              geom_line(aes(y = shape.97.5., col  = "97.5%")) +
              labs(title = "Shape parameter") +
              scale_colour_manual(name = "Quantiles", values = col.quantiles) +
              theme_piss()
            PissoortThesis::grid_arrange_legend(gloc.lin, gsc.lin, gsh.lin)
          }

          setProgress(1)


        })

        compare_methods_quantiles_plot(diff.boot.compp)
      }

    }

    else  print("Fill the box to allow for comparison") #ggplot()

  })


  output$datatable <- renderTable({

    time.fit <- fit()["time.fit"]
    time.boot <- boot()["time.boot"]

    df <- data.frame(Time.fit = time.fit, Time.boostrap = time.boot)
    colnames(df) <- c("Model fitting ",
                      "Bootstrapped intervals ")
    rownames(df) <- NULL
    #datatable(df, style = "bootstrap", selection = 'single')
    xtable::xtable(df)
  }, bordered = T, striped = T, hover = T)


  output$datatableFin <- DT::renderDataTable({

    df.boot <- boot()[["df.boot"]]

    df.q50 <- data.frame(Year = 1901:2016, Location = df.boot$loc.50.,
                         Scale = df.boot$scale.50., Shape = df.boot$shape.50.)
    #colnames(df.q50) = c("Year", "$\\mu \\ $", "$\\sigma \\quad$", "$\\xi \\quad$")

    datatable(round(df.q50,4), style = "bootstrap", selection = 'single',
              rownames = NULL, options = list(
                initComplete = JS(
                  "function(settings, json) {",
                  "$(this.api().table().header()).css({'background-color': '#33666C', 'color': '#fff'});",
                  # "$(this.api().table().first-child()).css({'background-color': '#33666C', 'color': '#fff'});",
                  "}" )))

  })


  'getPage' <- function(file = "information/infoNN.html") {
    return(includeHTML(file))
  }
  output$infoNN <- renderUI({ getPage() })



  #### 5th Application (Bayesian)


  data <- eventReactive(input$run, {
    # ## Start the progress bar.
    # withProgress(message = 'Gibbs sampling', value = 0, {

    # Create a Progress object
    progress <<- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Gibbs sampling", value = 0)


    data <- max_years$data

    fn <- function(par, data) -log_post1(par[1], par[2], par[3],
                                         par[4],rescale.time = T, data)
    param <- c(mean(max_years$df$Max), 0, log(sd(max_years$df$Max)), -0.1 )
    opt <- optim(param, fn, data = max_years$data,
                 method = "BFGS", hessian = T)
    opt

    Nbr.param <- length(opt$par)

    # Starting Values
    set.seed(input$seed)
    start <- list() ;   k <- 1  ;   n.chain <- input$start
    while(k < (n.chain+1)) { # starting values are randomly selected from a distribution
      # that is overdispersed relative to the target
      sv <- as.numeric(rmvnorm(1, opt$par, 50 * solve(opt$hessian)))
      svlp <- log_post1(sv[1], sv[2], sv[3], sv[4], max_years$data)
      if(is.finite(svlp)) {
        start[[k]] <- sv
        k <- k + 1
      }
    }
    mat_startvalues <- matrix(unlist(start), nrow = input$start, byrow = T)
    df_startvalues <- as.data.frame(mat_startvalues)

    set.seed(input$seed)
    iter.by.chain <- input$iterchain   ;  burnin = input$burnin

    # Handle the progress bar. See inside gibbs.trend.own()
    n.tot <- input$start * iter.by.chain
    'Progress.Shiny' <- function(detail = NULL) {
      progress$inc(amount = 1/n.tot, detail = detail)
    }

    # Handle the inputs for the Normal priors
    mu.mean.pr <- input$priormumean ;         mu.sd.pr <- input$priormusd
    mu1.mean.pr <- input$priormu1mean ;       mu1.sd.pr <- input$priormu1sd
    logsig.mean.pr <- input$priorlogsigmean ; logsig.sd.pr <- input$priorlogsigsd
    xi.mean.pr <- input$priorximean ;         xi.sd.pr <- input$priorxisd
    mean.vec <- c(mu.mean.pr, mu1.mean.pr,logsig.mean.pr, xi.mean.pr)
    sd.vec <- c(mu.sd.pr, mu1.sd.pr, logsig.sd.pr, xi.sd.pr)

    gibbs.trend <- PissoortThesis::gibbs.trend.own(start,
                                                   propsd = c(.5, 1.9, .15, .12),
                                                   iter = iter.by.chain,
                                                   burnin = burnin,
                                                   Progress.Shiny = Progress.Shiny, # Handles progress bar !
                                                   .mnpr = mean.vec, .sdpr = sd.vec)

    param.chain <- gibbs.trend$out.chain[, 1:Nbr.param]

    list(model = gibbs.trend, df_startvalues = df_startvalues, param.chain = param.chain)

  })

  output$plotPred1 <- renderPlot({
    #browser()

    mod <- data()[["model"]]
    from <-  input$from - 1900
    fut <-  input$fut
    by <- input$dens

    #browser()
    PissoortThesis::posterior_pred_ggplot(Data = max_years$df,
                                          Model_out.chain = mod$out.chain,
                                          from = from, x_coord = c(27, 35 + 0.02 * fut),
                                          n_future = fut, by = by)

  })


  output$plotPred2 <- renderPlot({

    mod <- data()[["model"]]
    from <-  input$from - 1900
    fut <-  input$fut
    by <- input$dens

    #browser()


    repl2 <- pred_post_samples(data = max_years$df, n_future = fut,
                               model_out.chain = mod$out.chain,
                               seed = input$seed, from = from)

    post.pred2 <- apply(repl2, 2,
                        function(x) quantile(x, probs = c(0.025,0.5,0.975)))
    hpd_pred <- as.data.frame(t(hdi(repl2)))


    if(fut == 0)  futur.dta <- NULL
    else   futur.dta <- repl2[sample(10, 1:nrow(repl2)), (ncol(repl2)-fut+1):ncol(repl2)]

    df.postpred2 <- data.frame(
      org.data = c(max_years$data[from:length(max_years$data)], futur.dta),
      q025 = post.pred2["2.5%",], q50 = post.pred2["50%",],
      q975 = post.pred2["97.5%",], year = input$from:(2016+fut),
      'data' = c(rep('original', length(max_years$data)-from+1), rep('new', fut)),
      hpd.low = hpd_pred$lower, hpd.up = hpd_pred$upper)

    col.interval <- c("2.5%-97.5%" = "red", "Median" = "blue2", "HPD 95%" = "green2",
                      "orange", "magenta")
    col.data <- c("original" = "cyan", "simulated" = "red", "orange", "magenta")

    g.ppd <- ggplot(df.postpred2) +
      geom_line(aes(x = year, y = q025, col = "2.5%-97.5%"), linetype = "dashed") +
      geom_line(aes(x = year, y = q50, col = "Median")) +
      geom_line(aes(x = year, y = q975, col =  "2.5%-97.5%"), linetype = "dashed") +
      geom_line(aes(x = year, y = hpd.low, col = "HPD 95%"), linetype = "dashed") +
      geom_line(aes(x = year, y = hpd.up , col =  "HPD 95%"), linetype = "dashed") +
      geom_vline(xintercept = 2016, linetype = "dashed", size = 0.4, col  = 1) +
      # scale_x_continuous(breaks = c(1900, 1950, 2000, 2016, 2050, 2100, 2131),
      #                    labels = c(1900, 1950, 2000, 2016, 2050, 2100, 2131) ) +
      scale_colour_manual(name = " PP intervals", values = col.interval) +
      geom_point(data = df.postpred2[1:116,],
                 aes(x = year, y = org.data), col = "black" ) +
      geom_point(data = df.postpred2[117:nrow(df.postpred2),],
                 aes(x = year, y = org.data), col = "orange" ) +
      scale_fill_discrete(name = "Data" ) + #, values = col.data) +
      labs(y = expression( Max~(T~degree*C)), x = "Year",
           title = "Posterior Predictive quantiles with observation + 116 years simulations") +
      theme_piss(size_p = 22, size_c = 19, size_l = 17,
                 theme = theme_minimal(),
                 legend.position =  c(0.91, 0.12))

    ## LEngth of the intervals
    length.quantil <- df.postpred2$q975 - df.postpred2$q025
    length.hpd <- df.postpred2$hpd.up - df.postpred2$hpd.low
    df.length.ci <- data.frame(quantiles = length.quantil,
                               hpd = length.hpd,
                               Year = df.postpred2$year)

    g.length <- ggplot(df.length.ci) +
      geom_line(aes(x = Year , y = quantiles), col = "red") +
      geom_line(aes(x = Year , y = hpd), col = "green2") +
      labs(title = "Intervals' lengths", y = "Length") +
      # scale_x_continuous(breaks = c(1900, 1950, 2000, 2050, 2100, 2131),
      #                    labels = c(1900, 1950, 2000, 2050, 2100, 2131) ) +
      geom_vline(xintercept = 2016, linetype = "dashed", size = 0.4, col  = 1) +
      theme(plot.title = element_text(size = 17, colour = "#33666C",
                                      face="bold", hjust = 0.5),
            axis.title = element_text(size = 10, colour = "#33666C", face="bold"))

    print(g.ppd)
    if (input$show){
      vp <- grid::viewport(width = 0.23,
                           height = 0.28,
                           x = 0.65,
                           y = 0.23)
      print(g.length, vp = vp)
    }



  })

  ## Gather the data for the traceplots
  traceplot.data <- reactive({
    mod <- data()[["model"]]

    chain.mix <- cbind.data.frame(mod$out.chain,
                                  iter.chain = rep( (input$burnin):(input$iterchain),
                                                    input$start))
    chain_mix_gg <- mixchains.Own(chain.mix, burnin = input$burnin)
    list(chain_mix_gg = chain_mix_gg)
  })

  ## Traceplots of the first parameters
  output$plot.chains <- renderPlot({
    chain_mix_gg <- traceplot.data()[["chain_mix_gg"]]

    title = "TracePlots of the generated Chains "
    grid_arrange_legend(chain_mix_gg$gmu, chain_mix_gg$gmutrend,
                        ncol = 2,
                        top = grid::textGrob(title,
                                             gp = grid::gpar(col = "#33666C",
                                                             fontsize = 25, font = 4)) )
  })

  ## Traceplots of the last parameters
  output$plot.chains2 <- renderPlot({
    chain_mix_gg <- traceplot.data()[["chain_mix_gg"]]

    grid.arrange(chain_mix_gg$glogsig, chain_mix_gg$gxi, ncol = 2)
  })

  ## Table of the starting values
  output$DTstart <- DT::renderDataTable({

    df <- data()[["df_startvalues"]]

    #rownames(df) <- replicate(1:input$start, paste("start ", i))
    colnames(df) <- c("mu", "mu1", "logsig", "xi")

    datatable(round(df,4), style = "bootstrap",
              selection = 'multiple', escape = F, rownames = NULL, options = list(
                initComplete = JS(
                  "function(settings, json) {",
                  "$(this.api().table().header()).css({'background-color': '#33666C', 'color': '#fff'});",
                  "}" ),
                dom = 't'))

  })

  output$accrates <- DT::renderDataTable({
    mod <- data()[["model"]]

    df_acc.rates <- matrix(unlist(mod$mean_acc.rates),
                           nrow = input$start, byrow = T) %>% t() %>%
      as.data.frame()
    mean.acc.rates <- colMeans(do.call(rbind, mod$mean_acc.rates))
    df <- cbind.data.frame(df_acc.rates, mean.acc.rates)

    colnames(df) <- c(paste0("start", 1:input$start), "Average")
    row.names(df) <- c("mu", "mu1", "logsig", "xi")

    datatable(round(df,4), style = "bootstrap",
              selection = 'multiple', escape = F, options = list(
                initComplete = JS(
                  "function(settings, json) {",
                  "$(this.api().table().header()).css({'background-color': '#33666C', 'color': '#fff'});",
                  "}" ),
                dom = 't')) %>%
      formatStyle( "Average",# target = 'row',
                   backgroundColor = "yellow"
      )
  })


  # Function to create mcmc.lists, useful for diagnostics on chains.
  'mc.listDiag' <- function(list, subset = c("mu0", "mu1", "logsig", "xi")) {
    if(length(list) <= 1 )
      resmc.list <- mcmc.list(mcmc(list[[1]][, subset]) )
    else {
      #browser()
      res <- list() # Initiialize list wherewe stock results
      res[[1]] <- mcmc(list[[1]][,subset])

      for (i in 2:length(list)) {

        res[[i]] <- mcmc( list[[i]][,subset] )
        # resmc.list <- mcmc.list(resmc.list,
        #                         res[[i]])
      }
      resmc.list <- mcmc.list(res)[,subset]
      #resmc.list <- lapply(res, mcmc.list )
    }
    return(resmc.list)
  }


  gg_gelman_reac <- reactive ({
    if(input$gelman) {

      mod <- data()[["model"]]

      gp.dat <- gelman.plot(mc.listDiag(mod$out.ind), autoburnin=F)
      df = data.frame(bind_rows(as.data.frame(gp.dat[["shrink"]][,,1]),
                                as.data.frame(gp.dat[["shrink"]][,,2])),
                      q = rep(dimnames(gp.dat[["shrink"]])[[3]],
                              each = nrow(gp.dat[["shrink"]][,,1])),
                      last.iter = rep(gp.dat[["last.iter"]], length(gp.dat)))
      df_gg <-melt(df, c("q","last.iter"), value.name = "shrink_factor")

      #browser()
      gg <-  ggplot(df_gg, aes(last.iter, shrink_factor, colour=q, linetype=q)) +
        geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
        geom_line() +
        geom_hline(yintercept = 1.1, colour = "green4", linetype = "dashed", size = 0.3) +
        scale_y_continuous(breaks = c(1, 1.1, 1.5, 2, 3, 4 ),
                           labels = c(1, 1.1, 1.5, 2, 3, 4 )) +
        #ggtitle("Gelman Rubin dignostic : R-hat Statistic") +
        facet_wrap(~variable,
                   labeller= labeller(.cols = function(x) gsub("V", "Chain ", x))) +
        labs(x="Last Iteration in Chain", y="Shrink Factor",
             colour="Quantile", linetype="Quantile",
             subtitle = "Gelman Rubin diagnostic : R-hat Statistic") +
        scale_linetype_manual(values=c(2,1)) +
        theme_piss() +
        theme(strip.text = element_text(size=15),
              plot.subtitle = element_text(size = 21, hjust = 0.5,
                                           colour = "#33666C", face = "bold"))
      return(gg)
    }
    else return(
      validate( need(input$gelman == T,
                     label = "Check the 'Gelman-R' box") )
    )
  })


  output$gelman <- renderPlot({ gg_gelman_reac() })


  gg_autocor <- reactive({

    if(input$autocor){
      param.chain <- data()[["param.chain"]]
      #browser()

      return(autocorr.plot(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")]  )) )
    }
    else return(
      validate( need(input$autocor == T,
                     label = "Check the 'autocorr' box") )
    )
  })
  output$autocorr <- renderPlot({ gg_autocor()   })


  gg_crosscor <- reactive({
    if(input$crosscor){
      param.chain <- data()[["param.chain"]]

      return(
        ggcorrplot(crosscorr(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")])),
                   hc.order = TRUE, type = "lower", lab = TRUE, title = "Cross-correlation",
                   ggtheme = PissoortThesis::theme_piss)
      )
    }
    else return(
      validate( need(input$crosscor == T,
                     label = "Check the 'Cross-corr' box") )
    )
  })
  output$crosscorr <- renderPlot({ gg_crosscor()  })



  geweke <- reactive({
    if(input$geweke){
      param.chain <- data()[["param.chain"]]

      return(
        geweke.plot(mcmc(param.chain), nbins = 20)
      )
    }
    else return(
      validate( need(input$geweke == T,
                     label = "Check the 'Geweke' box") )
    )
  })
  output$geweke <- renderPlot({ geweke()  })




  raftery <- reactive({
    if(input$raft){
      param.chain <- data()[["param.chain"]]

      return(
        raftery.diag(mcmc(param.chain[, c("mu0", "mu1", "logsig", "xi")]),
                     q=0.05, r=0.02, s=0.95)
      )
    }
    else return(
      validate( need(input$raft == T,
                     label = "Check the 'Raftery-Coda' box") )
    )
  })
  output$raft <- DT::renderDataTable({
    df <- as.data.frame(raftery()$resmatrix)

    datatable(round(df,4), style = "bootstrap",
              selection = 'multiple', escape = F, options = list(
                initComplete = JS(
                  "function(settings, json) {",
                  "$(this.api().table().header()).css({'background-color': '#33666C', 'color': '#fff'});",
                  "}" ),
                dom = 't'))
  })
  'getPage_raft' <- function(file = "information/info_raft.html") {
    return(includeHTML(file))
  }
  output$info_raft <- renderUI({ getPage_raft() })



  'getPage' <- function(file = "information/infobay.html") {
    return(includeHTML(file))
  }
  output$infobay <- renderUI({ getPage() })


}

shinyApp(ui, server)


#shiny::runApp(display.mode="showcase")

# options(shiny.trace = TRUE)
# options(shiny.fullstacktrace = TRUE)
# options(shiny.reactlog=TRUE)
