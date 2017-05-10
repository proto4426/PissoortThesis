library(shiny)
#options(shiny.error = browser)
library(PissoortThesis)
data("max_years")

ui <- fluidPage(

  titlePanel("Shiny app for master's thesis in extreme values"),
  selectInput("year",
              label = "Which year ? ",
              choices = as.character(seq(1901, 2016))),
  mainPanel(
    plotOutput("plot1", height = '600px', width = "900px")
  )
)
library(ggplot2)
library(plotly)
library(gridExtra)
library(grid)
library(highcharter)
library(dplyr)
library(tidyr)

load('C:\\Users\\Piss\\Documents\\LINUX\\Documents\\Extreme\\R resources\\IRM\\data.RData')





server <- function(input, output) {

  # hcbase <- reactive({
  #   # hcbase <- function() highchart()
  #   hc <- highchart()
  #
  #
  #   if (input$credits)
  #     hc <- hc %>% hc_credits(enabled = TRUE, text = "Highcharter", href = "http://jkunst.com/highcharter/")
  #
  #   if (input$exporting)
  #     hc <- hc %>% hc_exporting(enabled = TRUE)
  #
  #   if (input$theme != FALSE) {
  #     theme <- switch(input$theme,
  #                     null = hc_theme_null(),
  #                     economist = hc_theme_economist(),
  #                     dotabuff = hc_theme_db(),
  #                     darkunica = hc_theme_darkunica(),
  #                     gridlight = hc_theme_gridlight(),
  #                     sandsignika = hc_theme_sandsignika(),
  #                     fivethirtyeight = hc_theme_538(),
  #                     chalk = hc_theme_chalk(),
  #                     handdrwran = hc_theme_handdrawn()
  #     )
  #
  #     hc <- hc %>% hc_add_theme(theme)
  #   }
  #
  #   hc
  #
  # })

  output$plot1 <- renderHighchart({
    x <- TXTN_closed
     highchart() %>%
      hc_add_series(x[ x$year == input$year,], type = "scatter",
                    hcaes(x = Date, y = TX), name = 'TX' ) %>%
      hc_add_series(x[ x$year == input$year,], type = "scatter",
                    hcaes(x = Date, y = TN), name = 'TN' ) %>%
      hc_add_theme(hc_theme_flatdark()) %>%
     hc_title(text = "draggable points demo") %>%
     hc_xAxis(categories = month.abb) %>%
     hc_plotOptions(
       series = list(
         point = list(
           events = list(
             drop = JS("function(){
                    console.log(this.series)
                    window.data = _.map(this.series.data, function(e) { return e.y })
                    Shiny.onInputChange('inputname', data);
                    }"))
         )))

    # if(input$fit == 'loess') g <- g +
    #   stat_smooth(method = "loess", se = F, aes(colour = 'LOESS'))
  })

}
shinyApp(ui = ui, server = server)



