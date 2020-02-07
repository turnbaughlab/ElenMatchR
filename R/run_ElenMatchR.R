#' Run ElenMatchR Shiny app.
#'
#' Launches \code{ElenMatchR} in default browser.

run_ElenMatchR <- function() {
  
  appDir <- system.file("app",package = "ElenMatchR")
  if (appDir == "") {
    stop("Could not installation.", call = FALSE)
  }
  shiny::runApp(appDir)
}