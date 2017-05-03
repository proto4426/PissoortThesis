library(shiny)
#options(shiny.error = browser)


shinyUI(fluidPage(

  titlePanel("Shiny app for master's thesis in extreme values"),

  selectInput("max",
              "Do you want minima or maxima ? ", c('Max', 'Min') ),
  selectInput("fit",
              label = "Which fitting method ? ",
              choices = c("linear trend model" = "lm",
                          "nonparam trend model" =  "loess",
                          "broken linear trend" = "bl",
                          "broken linear + linear trend" = "blll",
                          "All 3 methods together" = 'all') ),

  mainPanel(
    plotOutput("plot1", height = '600px', width = "900px")
  )

))
