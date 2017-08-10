## app.R ##
library(shiny)


ui <- fluidPage(

  titlePanel("EVT thesis: Simulating GAM fit for the trend to visualize uncertainty"),

  fluidRow(
    column(3,
           selectInput("param",
                       label = "Which model ? ",
                       choices = c("Stationary" = "sta",
                                   "linear in location" =  "lin1",
                                   "linear in location an scale" = "lin2",
                                   "linear in location, scale and shape" = "lin3"),
           numericInput("bag", "Bootstrap resamples B ? ",
                        "1", min = "1", max = "500", step = 5 )
    ),


    fluidRow(
      column(3, offset = 1,
             numericInput("hidd", "Number of hidden layers ?",
                          "0", min = "0", max = "4" ),

             numericInput("beta", "Parameters for the beta prior ?" ,
                          "1", min = "2", max = "1000"  )
      ),
      fluidRow(
        column(3, offset = 1,
               selectInput("m",
                           label = "Which activation function ? ",
                           choices = c("Logistic sigmoid ?" = "logsig",
                                       "Hyperbolic tangent" =  "tanh")),

               numericInput("beta", "Parameters for the beta prior ?" ,
                            "1", min = "2", max = "1000"  )
        ) ),

      mainPanel(
        plotOutput("plot1", height = '500px', width = "750px")
      )
    ))
))


library(GEVcdn)
# library(foreach)
# library(doParallel)

library(PissoortThesis)
data("max_years")


x <- as.matrix(seq(1, length(max_years$data)))
y <- as.matrix(max_years$data)



server <- function(input, output) {


  if(input$bag == 1 ) {

  weights <-   PissoortThesis::gevcdn.fit2(x = x, y = y,
                                n.trials = 1,
                                n.hidden = input$hidd,
                                Th = ,
                                fixed = models[[i]]$fixed )

    }

  # cores <- detectCores()
  # cl <- makeCluster(cores[1]-1) #not to overload the computer, do not use all cores
  # registerDoParallel(cl)
  #
  # bag_par <- foreach(i = 1:M,
  #                    .packages = c("PissoortThesis", "GEVcdn"),
  #                    .verbose = T) %dopar% {
  #                      weights.on <- gevcdn.bag(x = x, y = y,
  #                                               iter.max = 100,
  #                                               fixed = c("shape"),
  #                                               iter.step = 10,
  #                                               n.bootstrap = 5,
  #                                               n.hidden = 2,
  #                                               sd.norm = .2)
  #                      parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)
  #
  #                      mean_of_list(parms.on)
  #                    }
}

shinyApp(ui, server)
