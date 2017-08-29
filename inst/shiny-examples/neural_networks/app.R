library(pander)
library(shiny)
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

library(PissoortThesis)
data("max_years")


#tanh activation function
'hyp.tan' <- function(x) ( exp(x) - exp(x) ) / ( exp(x) + exp(-x) )

# Create a function to easily compare all the methods. df contains the residuals for all quantiles
# See output$plot3
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




shinyApp(
ui <- tagList(
  navbarPage(  theme = shinytheme("united"),
    "EVT thesis of Antoine Pissoort",
    tabPanel("GEV-CDN : Neural Networks to fit flexible nonstationary Nonlinear GEV Models in parallel + prevent overfitting (bagging)",

  # Sidebar with sliders that demonstrate various available options
  sidebarPanel( #theme = shinythemes::shinytheme("united"),

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

    htmltools::p(actionButton("runfit","RUN GEV-CDN", icon("random") ),
      align = "center", width = 9 ),

    sliderInput("bag", "Bagging resamples B ? ",
                value = 1, min = 1, max = 500, step = 4 ),
    h5("1=no bagging. Only for nonstationary models only with hidden layer >0"),

    wellPanel(tags$h3("Bootstrap conf. intervals "),

      sliderInput("nboot", "Bootstrap resamples?",
                  value = 1, min = 1, max = 500, step = 4 ),
      selectInput("method",
                  label = "Bootstrap method ?",
                  choices = c("Residual" = "residual",
                              "Parametric" =  "parametric")),
      htmltools::p(actionButton("runboot","RUN Bootstrap intervals", icon("random") ),
        align = "center", width = 9 ),

    checkboxInput("comp", "Comparison of the resdiual and parametric methods ? ", FALSE)
    ),
      # numericInput("quant", "Quantiles ? ?" ,
      #              value = "6", min = "2", max = "1000"  ),
    wellPanel(tags$h5("Created by Antoine Pissoort"), tags$body("(", tags$a("Github",
                                    href="https://github.com/proto4426)"," | ",
                        tags$a("Linkedin", href="https://www.linkedin.com/in/antoine-pissoort-858b54113/"),")") )
    )
  ),

      mainPanel(
        tabsetPanel(#h3( "If bootstrap, look at console to follow computation's progress "),

        tabPanel(h5("Model fitting + parameters' intervals"),
                 # code("If bootstrap, look at console to follow computation's progress. \n
                 #       The computations are made in parralel using all your cores-1 by default.
                 #         But, computations with much bootstrap resamples could still take time."),
                 plotOutput("plot1", height = '500px', width = "800px"),
        br(),        br(),
        h4(strong("Computation times (sec.)")),
        tableOutput("datatable"),

        plotOutput("plot2", height = '700px', width = "600px")
        ),

        tabPanel("Interval's comparison + Summary Table",
                 # wellPanel(
                 #   h5( "Click on the  'comparison...' button below to see the graphs")
                 #   ),
                 h4(strong("Difference between residual and parametric Bootrstrap quantiles")),
                 plotOutput("plot3", height = '300px', width = "800px"),
                 br(),        br(),
                h5(strong("Parameters' Estimated values from the selected model")),
                DT::dataTableOutput("datatableFin", width = "800px")
        ),
        tabPanel("Informations", icon = icon("info-circle"),
                 htmlOutput("infoNN")  # See start of server()
        #tabPanel("Summary Table",  DT::dataTableOutput("datatableFin"))
      )
    )
  )
  )
  # footer = p(hr(), p("ShinyApp created by ", strong("{Antoine Pissoort}")),
  #                    p(("Code available on Github:"),a("https://github.com/proto4426",href="https://github.com/proto4426"))
  ) # navbarPage
),



server <- function(input, output) {


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


  # model <- list(Th = Th, beta.p = input$beta1, beta.q = input$beta2,
  #               n.hidden = hidd, fixed = )

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

 #browser()

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

      # # Create a Progress object
      # progress <<- shiny::Progress$new()
      # # Make sure it closes when we exit this reactive, even if there's an error
      # on.exit(progress$close())
      # progress$set(message = "Bagging", value = 0)


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
    # cl <- makeCluster(input$cores)
    #registerDoParallel(cl)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = M, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    t <- proc.time()

    incProgress(0.5)

    # 'Progress.Shiny' <- function(detail = NULL) {
    #   progress$inc(amount = 1/M, detail = detail)
    # }

    bag_par <- foreach(i = 1:M,
                       .packages = c("PissoortThesis", "GEVcdn"),
                       .options.snow = opts) %dopar% {

        weights.on <-  gevcdn.bag2(x = x, y = y,
                                   iter.max = 100,
                                                  fixed = fixed,
                                                  #iter.step = 10,
                                                  n.bootstrap = 2,
                                                  n.hidden = H,
                                                  sd.norm = sdnorm,
                                                  Th = Th,
                                                  silent = F)
                         parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)

                         #incProgress(1/M, detail = paste("Doing part", i))
                         # if (is.function(Progress.Shiny)) {
                         #   text <- paste0("Iteration ", i)
                         #   Progress.Shiny(#value = 1-(1/t*k),
                         #     detail = text)
                         # }
                         #progress$inc(amount = 1/M, detail = paste0("Iteration ", i))

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


  output$plot1 <- renderPlot({
    
    
    validate( need(input$runfit == T,
                   label = "Click on the 'RUN GEV-CDN' button ") )
    
    observeEvent(input$runfit == T , {
      df <- fit()[["df"]]
      
      output$plot1 <- renderPlot({
        
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
    })
    
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
  # cl <- makeCluster(input$cores)
  cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
  #registerDoParallel(cl)
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
                                                            #probs = c(0.1, 0.5, 0.9))
                                 list(loc = ci.lin$location.bootstrap,
                                      sc = ci.lin$scale.bootstrap,
                                      sh = ci.lin$shape.bootstrap
                                 )
                                # cat(i)
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

    list(#gg.loc.res = gg.boot.loc, gg.sc.res = gg.boot.loc, gg.sh.res = gg.boot.loc,
         df.boot = df.boot, time.boot = time.boot, meth = meth)

   })



  output$plot2 <- renderPlot({
    
    
    
    validate( need(input$runboot == T && input$runfit == T ,
                   label = "Click on the 'RUN GEV-CDN' and 'RUN Bootstrap' button ") )
    
    
    observeEvent(input$runboot == T, {
      meth <- boot()[["meth"]]
      df.boott <- boot()[["df.boot"]]   
      output$plot2 <- renderPlot({
        
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
                                    gp = gpar(fontsize = 24, font = 4, col ="black")))      })
    })
    
    
    
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
             PissoortThesis::grid_arrange_legend(gloc.lin, gsc.lin, gsh.lin )#,
                                                 # top ="Difference between residual and parametric Bootrstrap quantiles")
                                                 #                       gp = gpar(col ="orange", fontsize = 24, font = 2)) )
           }

           setProgress(1)

         })

         return(compare_methods_quantiles_plot(diff.boot.compp) )
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

}
)
#shinyApp(ui, server)
