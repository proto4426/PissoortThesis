# PissoortThesis
R package for thesis of Antoine Pissoort

## R code to install the package from GitHub
`devtools::install_github("proto4426/PissoortThesis", build_vignettes=T)`

Sometimes, you may have to use 

`devtools::install_github("proto4426/PissoortThesis", force=T)`

If you are not able to download the vignettes directly, please use the .tar.gz file containing the html's in the folder.

2. `library(PissoortThesis)`

## R code to install the package from a local repository:
1. `install.packages("path-to-PissoortThesis", repos = NULL, type="source")`
2. `library(PissoortThesis)`



# First visualisation : Shiny

After having loaded the package in your environement, you can run

`# Be sure to have grid, gridExtra, plotly and ggplot2 already installed`

`runExample()  # Present the yearly analysis visualizaiton for Block-Maxima (...)`


# Where can you find the DATA ? 

The project's data used for this thesis are confidential and have been distributed by the Institut Royal de Météorologie from Uccle. Hence, I cannot put it as public and you must ask me if you would like to handle the data. 

The `data` folder only contains the yearly data for use into the Shiny application. 
