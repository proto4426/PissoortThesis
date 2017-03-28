library(shiny)

load('/home/piss/Documents/Extreme/R resources/IRM/data1.Rdata')


ui <- fluidPage(

  titlePanel("Shiny app for master's thesis in extreme values"),

      selectInput("max",
                  "min or max ? ", c('Max', 'Min') ),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plot1", height = '600px', width = "900px")
    )

)

library(PissoortThesis)
library(ggplot2)
library(plotly)
library(gridExtra)
library(grid)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$plot1 <- renderPlot({
    x <- cbind.data.frame(max_years$df, Min = min_years$data)
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # ggplot(x) +
    #   geom_line(aes_string(x = 'Date', y = input$max))
     ggplot(x, aes_string(x = 'Year', y = input$max)) + geom_line() +
      geom_smooth(method='lm',formula=y~x, aes(colour = "Linear")) +
      # geom_line(data = max_years$df[max_years$df$Year %in% 1901:1975,],
      #           aes(x = Year, colour = "BrokenLinear",
      #               y = predict(lm(max_years$data[1:75] ~ max_years$df$Year[1:75]))),
      #           size = 1.5, linetype = "twodash") +
      # geom_line(data = max_years$df[max_years$df$Year %in% 1977:2016,],
      #           aes(x = Year, colour = "BrokenLinear",
      #               y = predict(lm(max_years$data[77:116] ~ max_years$df$Year[77:116]))),
      #           size = 1.5, linetype = "twodash") +
      stat_smooth(method = "loess", se = F, aes(colour = 'LOESS')) +
      labs(title = "Complete Serie of Annual TX in Uccle") +
      #theme_piss(20, 15) +
      theme(axis.line = element_line(color="#33666C", size = .45)) +
      scale_colour_manual(name="Trend",
                          values=c(Linear="blue", BrokenLinear="cyan", LOESS="red")) +
      #scale_linetype_manual( name = "Trend", values = c(BrokenLinear="twodash")
      theme(legend.position = c(.888, .152)) +
      theme(legend.title = element_text(colour="#33666C",
                                        size=12, face="bold"),
            legend.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = "white")) +
      guides(colour = guide_legend(override.aes = list(size = 2)))
    # print(ggplotly(g))
    # plot_ly(data = x, x = ~Date, y = ~input$stations, type = "scattergl") %>%
    #         #color = ~x, colors = c("green", "blue", "red"))
    #   layout(xaxis = xl, yaxis = yl, title = tit, legend = l)
    #
    # grid.arrange(g1, g2, ncol = 1,
    #              top = textGrob(expression(" Comparison of ns + solar cycle "),
    #                             gp = gpar(fontsize = 22, font = 3,
    #                                       col ="#33666C")))
  })

}



shinyApp(ui = ui, server = server)


#runApp('/home/piss/Documents/Shiny/appTest')

