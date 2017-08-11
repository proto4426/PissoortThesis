load("/home/proto4426/Documents/Thesis/Extreme/R resources/IRM/data1.Rdata")

library(tidyverse)
library(PissoortThesis)
library(grid)

##### Draw the plot of the beggining of Chapter 2


p <- ggplot(TXTN_closed) + geom_density(aes(x = TX), fill = "grey", alpha = .5) +
  labs(title = "Distribution of maximum daily temperatures",
       x = expression( bold(Maximum~T~degree*C)))
d <- ggplot_build(p)$data[[1]]

p <- p + geom_area(data = subset(d, x > 30), aes(x=x, y=y), fill="red", alpha = .8) +
  geom_segment(x=30, xend = 30, y=0, yend=approx(x = d$x, y = d$y, xout = 30)$y,
               colour="blue", size=1.2) +
  theme_piss(theme = theme_classic()) +
  coord_cartesian(xlim = c(-7, 39)) +
  geom_segment(aes(x = 32.5, xend = 32.5, y = .013, yend = .00295),
               arrow = arrow(length = unit(0.5, "cm")), col = "red") +
  annotate("text", label = "u",
           x = 30, y = .007, col = "blue", size = 7, fontface = 3)
p


gg_above30 <- ggplot(data=above_30, aes(TX),) +
  geom_density(alpha = .3, size=1, fill = "red", col = "red") +
  labs(title = expression(bold( Above~30~degree*C))) +
  theme_piss(theme = theme_classic()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color="#33666C", size = .5),
        axis.text.x = ) +
  theme(axis.title = element_blank())



vp <- viewport(width = 0.26,
               height = 0.43,
               x = 0.86,
               y = 0.57)
p
print(gg_above30, vp = vp)
