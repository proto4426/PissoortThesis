library(shiny)
#options(shiny.error = browser)


shinyUI(fluidPage(

  titlePanel("EVT thesis: Simulating GAM fit for trend to visualize uncertainty"),

  fluidRow(
    column(3,
  numericInput("level", withMathJax(helpText("Which level $\alpha$ ? (in %) ")),
               "5", min = "0", max = "100" ),

  numericInput("sim", "Howmuch Simulations M ? ",
               "50", min = "2", max = "1000" )
   ),  width = "100px",


  fluidRow(
    column(3, offset = 1,
  numericInput("seed", "Set the seed ",
               "1234", min = "1", max = "1000000000" ),

  numericInput("draws", "Howmuch simulations to draw ? (must be < M)",
               "50", min = "2", max = "1000"  )
    ),

  mainPanel(
    plotOutput("plot1", height = '600px', width = "900px")
  )
))
))
