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

### 4. Neural Networks (GEV-CDN)
Present the Conditional Density Netowrks applied to the GEV (see Sections **3.4** and  **6.3**). This app provides convenient visualization to see the effect of all the (hyper)parameters of the model(s).  
```coffee
runExample('neural_networks') 
```
Parralel computing is used for efficiency. Click on the "informations" tab in the application for more explanations. Here is a very quick summary...
![gap_test](https://github.com/proto4426/LaTeX_new/blob/master/gif/NN_small.gif)


### 5. Bayesian Analysis
Present the Bayesian Analysis applied in Chapter **7** in the thesis. This app also provides convenient tools to visualize the effects of changing the (hyper)parameters of the model(s), diagnostics, Predictive Posterior, ... 
```coffee
runExample('Bayesian') 
```
Click on the "informations" tab in the application for more explanations. Here is a very quick GIF summary :
![gap_test](https://github.com/proto4426/LaTeX_new/blob/master/gif/Bayes.gif)
**UPDATE : ** C++ backend for Gibbs sampling is now enabled, relying on the `Rcpp` R package. An  $\approx$**40 times increase of computational time efficiency** is observed. Comparisons can be found inside the application for the generated chains.


### All together in a dashboard 
All applications created above in a smooth dashboard. Allows to put quite more tools. 
```coffee
runExample('All_dashboard') 
```
can be directly accessed by clicking on this [link](https://proto4426.shinyapps.io/All_dashboard/).
However, it is recommended to visualize the *Bayesian* and *Neural Networks* applications individually since the dashboard's layout is not optimal with the panel divided in tabsets.



# Where can you find the DATA ? 

The project's data used for this thesis are confidential and have been provided by the Institut Royal de Météorologie from Uccle. Hence, I cannot put it as public and you must ask me if you want to obtain the data. 

The `data` folder only contains the yearly data to allow use into the Shiny application. 
