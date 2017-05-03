library(shiny)
#options(shiny.error = browser)


shinyUI(fluidPage(

  titlePanel("Thesis in EVT : Simulating a GAM trend to visualize uncertainty"),

  numericInput("sim", "Howmuch Simulations M to make ? ",
               "50", min = "2", max = "1000" ),

  numericInput("seed", "Set the seed ",
               "1234", min = "1", max = "1000000000" ),

  numericInput("draws", "Howmuch simulations to draw ? (must be < M)",
               "50", min = "2", max = "1000" , ),


  mainPanel(
    plotOutput("plot1", height = '600px', width = "900px")
  )

))
