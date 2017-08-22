library(coda)
library(htmltools)
library(shiny)
library(shinythemes)
library(shinydashboard)

library(gridExtra)
library(grid)
library(DT)
library(tidyverse)
library(reshape2)

library(mvtnorm)
library(HDInterval)
library(ggjoy)
library(viridis)
library(plotly)
library(ggcorrplot)


library(PissoortThesis)
data("max_years")


shinyApp(
 ui <- tagList(
 navbarPage(
   theme = shinytheme("spacelab"),
   "EVT thesis of Antoine Pissoort",
    tabPanel("Bayesian Analysis",
     sidebarPanel(
       wellPanel(tags$h2("Model"),
       sliderInput("iterchain", "Number of iterations by chains ",
                   value = 500, min = 10, max = 1e4, step = 20),
       numericInput("seed", "Set the seed ",
                   value = 123, min = 1, max = 1e15)
       ),
       wellPanel(tags$h2("Priors"),
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
       tags$h5("Take them uninformative by default ")
       ),

       htmltools::p(actionButton("run","RUN Gibbs Sampler", icon("random") ),
         align = "center", width = 9 ),

       wellPanel(tags$h2("Diagnostics"),
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
              )
       ),
       br(),
       wellPanel(tags$h2("Posterior Predictive"),
       sliderInput("from", "Start at year ",
                value = 1901, min = 1901, max = 2016, step = 1 ),
       sliderInput("fut", "Prediction in the future ? ",
                value = 1, min = 0, max = 250, step = 5 ),
       numericInput("dens", "Show densities every which years ?",
                value = "10", min = "1", max = "30" ),
       #h4("years"),
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
          plotOutput("plot1", height = '800px', width = "800px"),
          br(),        br(),
          plotOutput("plot2", height = '500px', width = "800px"),

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
                # box(title = "Showed MCMC Diagnostics",
                #     status = "warning", solidHeader = TRUE, # collapsible = TRUE,
                #     background = "light-blue",
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
                  )   ## ============================================================
                     # fluidRow(
                    #   column(10, plotOutput(outputId="autocorr", width="450px",height="400px")),
                    #   column(4, plotOutput(outputId="crosscorr", width="250px",height="400px"))
                    #)
                  )
                ),
                tabPanel("Informations", icon = icon("info-circle"),
                         htmlOutput("info")#,
                         #actionLink("link_to_info", "Link Informations")
                )
            )  # tabsetPanel
        )  # mainPanel
       ), # tabPanel
   br(), br(), br(), br(),
   footer = htmltools::p( hr(), htmltools::p("ShinyApp created by ", strong("Antoine Pissoort"), "(",
                                  a(icon("linkedin"), "Linkedin",
                                    href = "https://www.linkedin.com/in/antoine-pissoort-858b54113/"), ")", align="center",width=4),
                          htmltools::p(("Code available on "), a(icon("github"),"Github",
                                          href="https://github.com/proto4426"),align="center",width=4)
           )
    )#, # navbarPage
#  tags$footer(title="Your footer here", align = "right", style = "
# position:absolute;
# bottom:0;
# width:100%;
# height:50px; /* Height of the footer */
# color: white;
# padding: 10px;
# background-color: black;
# z-index: 1000;"
#  )

 ), # tagList


server <- function(input, output) {

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

   output$plot1 <- renderPlot({
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

   output$plot2 <- renderPlot({

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


   'getPage' <- function(file = "information/info.html") {
     return(includeHTML(file))
   }
   output$info <- renderUI({ getPage() })

   ## Link to info ? https://stackoverflow.com/questions/34315485/linking-to-a-tab-or-panel-of-a-shiny-app
   # observeEvent(input$link_to_info, {
   #   newvalue <- "B"
   #   updateTabItems(session, "panels", newvalue)
   # })


})

# shiny::runApp(display.mode="showcase")
#
# options(shiny.trace = TRUE)
# options(shiny.fullstacktrace = TRUE)
# options(shiny.reactlog=TRUE)

