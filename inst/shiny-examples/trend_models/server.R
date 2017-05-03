# load('C:\\Users\\Piss\\Documents\\LINUX\\Documents\\Extreme\\R resources\\IRM\\min_years.Rdata')
# save(min_years, file = "C:\\Users\\Piss\\Documents\\LINUX\\PissoortRepo\\PissoortThesis\\data\\min_years.RData")
# load("C:\\Users\\Piss\\Documents\\LINUX\\PissoortRepo\\PissoortThesis\\data\\data_shiny_yearly.RData")

data("min_years")
data("max_years")

library(PissoortThesis)
library(ggplot2)
library(plotly)
library(gridExtra)
library(grid)



shinyServer(function(input, output) {

  output$plot1 <- renderPlot({
    x <- cbind.data.frame(max_years$df, Min = min_years$data)


    g <- ggplot(x, aes_string(x = 'Year', y = input$max)) + theme_bw() + geom_line() + geom_point() +
      theme(plot.title = element_text(size = 22, hjust=0.5,
                                      colour = "#33666C", face="bold"),
            axis.title = element_text(face = "bold", size= 18,
                                      colour = "#33666C"),
            axis.line = element_line(color="#33666C", size = .45),
            legend.position = c(.888, .152),
            legend.background = element_rect(colour = "black"),
            legend.key = element_rect(fill = "white") ) +
      scale_colour_manual(name="Trend",
                          values=c(Linear="blue", BrokenLinear="cyan", LOESS="red")) +
      guides(colour = guide_legend(override.aes = list(size = 4)))


    g_linear <- geom_smooth(method='lm',formula=y~x, aes(colour = "Linear"))

    g_smooth <-  stat_smooth(method = "loess", se = F, aes(colour = 'LOESS'))

    g_bl_max <-  list(geom_line(data = max_years$df[max_years$df$Year %in% 1901:1975,],
                                aes(x = Year, colour = "BrokenLinear",
                                    y = predict(lm(max_years$data[1:75] ~ max_years$df$Year[1:75]))),
                                size = 1.5, linetype = "twodash"),
                      geom_line(data = max_years$df[max_years$df$Year %in% 1977:2016,],
                                aes(x = Year, colour = "BrokenLinear",
                                    y = predict(lm(max_years$data[77:116] ~ max_years$df$Year[77:116]))),
                                size = 1.5, linetype = "twodash")     )

    g_bl_min <- list(geom_line(data = min_years$df[min_years$df$Year %in% 1901:1975,],
                               aes(x = Year, colour = "BrokenLinear",
                                   y = predict(lm(min_years$data[1:75] ~ min_years$df$Year[1:75]))),
                               size = 1.5, linetype = "twodash"),
                     geom_line(data = min_years$df[min_years$df$Year %in% 1977:2016,],
                               aes(x = Year, colour = "BrokenLinear",
                                   y = predict(lm(min_years$data[77:116] ~ min_years$df$Year[77:116]))),
                               size = 1.5, linetype = "twodash")     )


    if(input$fit == 'lm')  g <- g + g_linear
    if(input$fit == 'loess') g <- g + g_smooth


    ## Broken linear trends

    if (input$max == "Max" & input$fit == "bl")
      g <- g + g_bl_max

    if (input$max == "Min" & input$fit == "bl")
      g <- g + g_bl_min


    ## Broken linear and linear trend

    if (input$max == "Max" & input$fit == "blll")
      g <- g + g_linear + g_bl_max

    if (input$max == "Min" & input$fit == "blll")
      g <- g + g_linear + g_bl_min


    ### ALL the methods

    if (input$max == "Max" & input$fit == "all")
      g <- g + g_linear + g_bl_max + g_smooth

    if (input$max == "Min" & input$fit == "all")
      g <- g + g_linear +  g_bl_min + g_smooth



    if(input$max == "Min") g <- g + labs(title = "Complete Serie of Annual TN in Uccle")
    else g <- g + labs(title = "Complete Serie of Annual TX in Uccle")
    g
  })

})
#shinyApp(ui = ui, server = server)
