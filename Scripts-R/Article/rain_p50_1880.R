library(data.table)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

library(PissoortThesis)

# Apply the created theme to all ggplots without having to specify it !
theme_set(PissoortThesis::theme_piss())

######### Read Data ###########
###############################

# From 1833. Free dataset from the KNMI
rain_1880 <- read.csv("/home/proto4426/Documents/master_Thesis/Extreme/R resources/IRM/P50_Uccle_1880.csv",
                      sep = " ")
names(rain_1880)[1] <- "Date"

rain_1880$day <- substr(rain_1880$Date,7,8)
rain_1880$month <- as.numeric(substr(rain_1880$Date,5,6))
rain_1880$year <- substr(rain_1880$Date,1,4)


# Retrieve seasons with our function. Based on meteorological seasons
rain_1880$season <- sapply(rain_1880$month,
                             function(x) PissoortThesis::func_season(x))

rain_1901 <- rain_1880[rain_1880$year > 1900,]


ggplot(rain_1901) + geom_histogram(aes(x = RR))



ggplot(data=rain_1901, aes(group = month)) +
  geom_boxplot(aes(x = month, y = RR)) +
  theme_bw()

## Violin-plots
dodge <- position_dodge(width = 0.4)
gv1 <- ggplot(rain_1901, aes(x = season, y = RR)) +
  geom_jitter(color='red', size = .6, alpha=0.99,width = 0.2) +
  geom_violin(fill = "lightseagreen", alpha=0.7, draw_quantiles = T,
              position = dodge, width = 1.8) +
  geom_boxplot(width=.06, outlier.colour=NA, position = dodge) +
  labs(title = 'Violin-plots for  by seasons',
       x = "Season", y = "mm") +
  theme_piss( size_p = 16)

## Density plots
ggplot(data = rain_1901, aes(RR, colour = as.factor(month))) +
  geom_density(size = 1.1) +
  scale_color_discrete() +
  theme_bw()
# !! same smoothing factor for all densities

summer <- rain_1901[rain_1901$season == "Summer", ]
spring <- rain_1901[rain_1901$season == "Spring", ]
winter <- rain_1901[rain_1901$season == "Winter", ]
autumn <- rain_1901[rain_1901$season == "Autumn", ]

m_summer <- mean(summer$RR)
m_spring =  mean(spring$RR)
m_winter <- mean(winter$RR)
m_autum <- mean(autumn$RR)

gd1 <- ggplot(data=rain_1901, aes(RR, fill = season, colour=season)) +
  geom_density(alpha = .1, size=1.1) +
  scale_fill_brewer(palette = "Set1" )+
  scale_color_brewer(palette= "Set1") +
  geom_hline(yintercept=0, colour="white", size=1.1) +
  labs(title = 'Kernel Densities for daily Max. t°c by seasons',
       y = "Density", x = expression( Maximum~T~degree*C)) +
  theme_piss(legend.position = c(0.9, .82), size_p = 16) +
  geom_vline(xintercept = m_summer, colour = "darkgreen", linetype = 2) +
  geom_vline(xintercept = m_spring, colour = "blue", linetype = 2) +
  geom_vline(xintercept = m_winter, colour = "violet", linetype = 2) +
  geom_vline(xintercept = m_autum, colour = "red", linetype = 2)


## violin and density plots together
grid.arrange(gv1, gd1, nrow = 1)


###########################


# block length : the usual method is one 1 year
list_rain_by_years <- split(rain_1901, rain_1901$year)
# Then, we have 116 years of data ! (1901 to 2016)


## Retrieve the max (TN) in each year
max_years <- PissoortThesis::yearly.extrm(list.y = list_rain_by_years, tmp = "RR")


## Plot the yearly maxima together with some "usual" fitting methods :
#linear regression, LOESS and broken linear regression
lm1 <- lm(max_years$data ~ max_years$df$Year)
lm1_1 <- lm1$coefficients[1]
lm1_2 <- lm1$coefficients[2]
summary(lm1)

gg_trends <- ggplot(data = max_years$df,aes(x=Year,y=Max)) +
  geom_line() + geom_point() +
  geom_smooth(method='lm',formula=y~x, aes(colour = "Linear")) +
  stat_smooth(method = "loess", se = F, aes(colour = 'LOESS')) +
  labs(title = "Complete Serie of Annual TX in Uccle",
       y = expression( Max~(T~degree*C))) +
  theme(axis.line = element_line(color="#33666C", size = .45)) +
  scale_colour_manual(name="Trend", values=c(Linear="blue",
                                             BrokenLinear="cyan",
                                             LOESS="red")) +
  theme_piss(legend.position = c(.92, .12), size_p = 23) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
gg_trends


# Histogram and density plots of the annual maxima
ggplot(data = max_years$df, aes(x=Max)) +
  geom_histogram() +
  theme_minimal()
ggplot(data = max_years$df, aes(x=Max)) +
  geom_density() +
  theme_bw()
