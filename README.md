# PissoortThesis
R package for the master's thesis in statistics of Antoine Pissoort at Université Catholique de Louvain

## R code to install the package from GitHub
Install the `devtools` package if you don't have it yet. Then 

```coffee
devtools::install_github("proto4426/PissoortThesis", build_vignettes=T)
```

Sometimes, you may have to use 

```coffee
devtools::install_github("proto4426/PissoortThesis", force=T)
```

If you are not able to download the vignettes directly, please use the .tar.gz file containing the html's in the folder.

```coffee
library(PissoortThesis)
```

## R code to install the package from a local repository:

```coffee
install.packages("path-to-PissoortThesis", repos = NULL, type="source")
library(PissoortThesis)
```


# First visualisation : Shiny

After having loaded the package in your environement, with one signle line of code you can run the shiny applications in your local environment : 

### 1.) GEV Distribution 
Present the GEV distribution and the dependence with its parameters
```coffee
runExample('GEV_distributions') 
# Make sure to have grid, gridExtra, plotly and ggplot2 already installed
```
![gap_test](https://github.com/proto4426/LaTeX_new/blob/master/gif/gev_distrib.gif)
*Visual "problems" with the EV Weibull upper end point (red) is due to a `ggplot2` mispecification in y-scale*


### 2. Models for the Trend
Present the yearly analysis visualizaiton for yearly Maxima (see Section 5.2.2)
```coffee
runExample('trend_models')  
```
![gap_test](https://github.com/proto4426/LaTeX_new/blob/master/gif/trend_models.gif)

### 3. Splines draws with GAM 
Present the simulation study of the GAM model with splines (see Section 5.2.3)
```coffee
runExample('splines_draws') 
```
![gap_test](https://github.com/proto4426/LaTeX_new/blob/master/gif/splines.gif)


### 4. All together in a dashboard 
All applications created above in a smooth dashboard. Allows to put quite more tools. 
```coffee
runExample('All_dashboard') 
```
can be directly accessed by clicking on this [https://proto4426.shinyapps.io/All_dashboard/](link)



# Where can you find the DATA ? 

The project's data used for this thesis are confidential and have been distributed by the Institut Royal de Météorologie from Uccle. Hence, I cannot put it as public and you must ask me if you would like to handle the data. 

The `data` folder only contains the yearly data to allow use into the Shiny application. 
