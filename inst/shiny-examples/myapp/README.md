## To be able to load the applications, you must have the `.RData` file...

Please contact the maintainer to do so, as the data are `confidential´.



# How to Use

## The easiest way (recommended)

After you have installed and loaded the package, simply run the command
```{r}
PissoortThesis::runExample()
```
and the application will appear in Rstudio.

## Or you can also...

Many ways other ways to download and run it:
  ```R
library(shiny)

# Easiest way is to use runGitHub
runGitHub("Shiny_visual_app", "rstudio")

# Run a tar or zip file directly
runUrl("https://github.com/proto4426/PissoortThesis/tree/master/Shiny_visual_app/archive/master.tar.gz")
runUrl("https://github.com/proto4426/PissoortThesis/tree/master/Shiny_visual_app/archive/master.zip")
```

Or you can clone the git repository, then use `runApp()`:

  ```R
# First clone the repository with git. If you have cloned it into
# ~/shiny_example, first go to that directory, then use runApp().
setwd("~/Shiny_visual_app")
runApp()
```


To run a Shiny app from a subdirectory in the repo or zip file, you can use the `subdir` argument. This repository happens to contain another copy of the app in `inst/shinyapp/`.

```R
runGitHub("Shiny_visual_app", "rstudio", subdir = "inst/shinyapp/")

runUrl("https://github.com/proto4426/PissoortThesis/tree/master/Shiny_visual_app/archive/master.tar.gz",
       subdir = "inst/shinyapp/")
```
