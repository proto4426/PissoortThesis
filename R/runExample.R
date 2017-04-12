#' @export runExample
#' @title Shiny application
#' @examples
#' ValUSunSSN::runExample() # and you should see the app
'runExample' <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "PissoortThesis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `ValUSunSSN`.", call. = F)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
