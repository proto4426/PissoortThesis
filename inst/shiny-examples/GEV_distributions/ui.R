library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Shiny app for Val-U-Sun project! Compare 2 stations"),


  # Sidebar with a slider input for the number of bins
  fluidRow(
    column(3,
           numericInput("ksi1", "Which value of shape parameter for Weibull-type",
                        "-0.5", min = "-100", max = "0" ),
           #numericInput("ksi2", "Which value of shape for 2nd density","0", min = "-100", max = "100" ),
           numericInput("ksi3", "Which value of shape parameter for Fréchet-type",
                        "0.5", min = "0", max = "100" )
    ),
    width = "100px",

    column(3, offset = 1,
           numericInput("mu", "Which location parameter ?",
                        "0", min = "-100000", max = "100000" ),
           numericInput("sig", "Which scale parameter ?",
                        "1", min = "0", max = "10000" )
    )
           ,

    # Show a plot of the generated distribution
    #mainPanel(
    plotOutput("plot1", height = '650px', width = "900px")
    #)

  )
))
