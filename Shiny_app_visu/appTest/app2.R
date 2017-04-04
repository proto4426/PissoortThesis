library(shiny)
#options(shiny.error = browser)
load('/home/piss/Documents/Extreme/R resources/IRM/data1.Rdata')


ui <- fluidPage(

  titlePanel("Shiny app for master's thesis in extreme values"),
  selectInput("year",
              label = "Which year ? ",
              choices = as.character(seq(1901, 2016))),
  mainPanel(
    plotOutput("plot1", height = '600px', width = "900px")
  )
)
library(PissoortThesis)
library(ggplot2)
library(plotly)
library(gridExtra)
library(grid)
library(highcharter)
library(dplyr)
library(tidyr)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$plot1 <- renderHighchart({
    x <- TXTN_closed
   hc <-  highchart() %>%
      hc_add_series(x[ x$year == input$year,], type = "scatter",
                    hcaes(x = Date, y = TX), name = 'TX' ) %>%
      hc_add_series(x[ x$year == input$year,], type = "scatter",
                    hcaes(x = Date, y = TN), name = 'TN' ) %>%
      hc_add_theme(hc_theme_flatdark())
  hc
    # if(input$fit == 'loess') g <- g +
    #   stat_smooth(method = "loess", se = F, aes(colour = 'LOESS'))
  })

}
shinyApp(ui = ui, server = server)



