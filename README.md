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
# Be sure to have grid, gridExtra, plotly and ggplot2 already installed
```
Example of what this application does : 
![gap_test](https://github.com/proto4426/Thesis/blob/master/gapminder2.gif)


### 2. Models for the Trend
Present the yearly analysis visualizaiton for Block-Maxima (see Section 6.2.2)
```
runExample('trend_models')  
```

### 3. Splines draws with GAM 
Present the simulation study of the GAM model with splines (see Section 6.2.3)
```
runExample('splines_draws') 
```


# Where can you find the DATA ? 

The project's data used for this thesis are confidential and have been distributed by the Institut Royal de Météorologie from Uccle. Hence, I cannot put it as public and you must ask me if you would like to handle the data. 

The `data` folder only contains the yearly data for use into the Shiny application. 
