#' @export runExample
#' @title Shiny application
#' @param example  chooses the example you want. Run runExample() to have the list
#' @examples
#' library(PissoortThesis)
#' runExample('Name of the app') # and you will see the app you want
#'# You juste have to run runExample('') with whatever incorrect name to display the names
#'# of the existing apps
#'
'runExample' <- function(example) {
  # locate all the shiny app examples that exist
  validExamples <- list.files(system.file("shiny-examples", package = "PissoortThesis"))

  validExamplesMsg <-
    paste0(
      "Valid examples are: '",
      paste(validExamples, collapse = "', '"),
      "'")

  # if an invalid example is given, throw an error
  if (missing(example) || !nzchar(example) ||
      !example %in% validExamples) {
    stop(
      'Please run `runExample()` with a valid example app as an argument.\n',
      validExamplesMsg,
      call. = FALSE)
  }

  # find and launch the app
  appDir <- system.file("shiny-examples", example, package = "PissoortThesis")
  shiny::runApp(appDir, display.mode = "normal")
}

# 'runExample' <- function() {
#   appDir <- system.file("shiny-examples", "myapp", package = "ValUSunSSN")
#   if (appDir == "") {
#     stop("Could not find example directory. Try re-installing `ValUSunSSN`.", call. = F)
#   }
#
#   shiny::runApp(appDir, display.mode = "normal")
# }
