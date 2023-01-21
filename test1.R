foo <- function() {
  message("one")
  Sys.sleep(0.5)
  message("two")
}

runApp(shinyApp(
  ui = fluidPage(
    shinyjs::useShinyjs(),
    actionButton("btn","Click me"),
    textOutput("text")
  ),
  server = function(input,output, session) {
    observeEvent(input$btn, {
      withCallingHandlers({
        shinyjs::html("text", "")
        foo()
      },
      message = function(m) {
        shinyjs::html(id = "text", html = m$message, add = TRUE)
      })
    })
  }
))