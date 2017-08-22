library(shiny)
library(shinydashboard)

header <- dashboardHeader( title = "Extreme Values' Analysis in Uccle")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("GEV distribution", tabName = "gev", icon = icon("dashboard"),
             badgeLabel = "chap.1", badgeColor = "green"),
    menuItem("Trend Models", icon = icon("thermometer-1"), tabName = "trend",
             badgeLabel = "chap.5", badgeColor = "green"),
    menuItem("Splines GAM", icon = icon("coffee"), tabName = "splines",
             badgeLabel = "chap.5", badgeColor = "red"),
    menuItem("Neural Networks", icon = icon("plane"), tabName = "nn",
             badgeLabel = "chap.6", badgeColor = "red"),
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
  tabItem(tabName = "nn",
          navbarPage(
            theme = shinytheme("united"),
            "EVT thesis of Antoine Pissoort",
            tabPanel("GEV-CDN : Neural Networks to fit flexible nonstationary Nonlinear GEV Models in parallel + prevent overfitting (bagging)",
                     sidebarPanel(
                       wellPanel(tags$h3("Hyperparameters "),
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
                       p(actionButton("runfit","RUN GEV-CDN", icon("random") ),
                         align = "center", width = 9 ),
                       sliderInput("bag", "Bagging resamples B ? ",
                                   value = 1, min = 1, max = 500, step = 4 ),
                       h5("1=no bagging. Only for nonstationary models with hidden layer >0"),
                       wellPanel(tags$h3("Bootstrap conf. intervals "),
                                 sliderInput("nboot", "Bootstrap resamples?",
                                             value = 1, min = 1, max = 500, step = 4 ),
                                 selectInput("method",
                                             label = "Bootstrap method ?",
                                             choices = c("Residual" = "residual",
                                                         "Parametric" =  "parametric")),
                                 p(actionButton("runboot","RUN Bootstrap intervals", icon("random") ),
                                   align = "center", width = 9 ),
                                 checkboxInput("comp", "Comparison of the resdiual and parametric methods ? ", FALSE)
                       ),
                       wellPanel(tags$h5("Created by Antoine Pissoort"), tags$body("(", tags$a("Github",
                                  href="https://github.com/proto4426)"," | ",
                                  tags$a("Linkedin", href="https://www.linkedin.com/in/antoine-pissoort-858b54113/"),")") )
                       )
                     ),
                     mainPanel(
                       tabsetPanel(#h3( "If bootstrap, look at console to follow computation's progress "),
                         tabPanel(h5("Model fitting + parameters' intervals"),
                                  plotOutput("plotFitNN", height = '500px', width = "800px"),
                                  br(),        br(),
                                  h4(strong("Computation times (sec.)")),
                                  tableOutput("datatable"),

                                  plotOutput("plotBoot", height = '700px', width = "600px")
                         ),
                         tabPanel("Interval's comparison + Summary Table",
                                  plotOutput("plot3", height = '300px', width = "780px"),
                                  br(),        br(),
                                  h5(strong("Parameters' Estimated values from the selected model")),
                                  DT::dataTableOutput("datatableFin")
                         ),
                         tabPanel("Informations", icon = icon("info-circle"),
                                  htmlOutput("infoNN")  # See start of server()
                                  #tabPanel("Summary Table",  DT::dataTableOutput("datatableFin"))
                         )))
            )
          ) # navbarPage, h2("Refer to Section ", strong("3.4"), "and", strong("6.3"), "for more information")
  )
 )
)


ui <- dashboardPage(header, sidebar, body)


library(evd)
library(PissoortThesis)
library(plotly)
library(gridExtra)
library(grid)
library(htmltools)
library(shinythemes)
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


  output$plotFitNN <- renderPlot({

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



  output$plotBoot <- renderPlot({

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


  output$plot3 <- renderPlot({

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


  'getPage' <- function(file = "~/inst/shiny-examples/neural_networks/information/infoNN.html") {
    return(includeHTML(file))
  }
  output$infoNN <- renderUI({ getPage() })

}

shinyApp(ui, server)



