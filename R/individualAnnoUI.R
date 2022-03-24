DEAnnotationReportServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {

      output$DE_AnnotationReport <- renderDataTable({

        dbareport()

      })

    }
  )
}

DEAnnotationReportUI <- function(id) {
  ns <- NS(id)
  tagList(
    dataTableOutput(ns("DE_AnnotationReport"))
  )
}
