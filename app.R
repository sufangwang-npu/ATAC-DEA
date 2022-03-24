options(shiny.maxRequestSize = 30*1024^2)

source("server.R")
source("ui.R")

options(shiny.sanitize.errors = FALSE)
setwd("/srv/shiny-server/ATAC-DEA")

annoData_Hsapines <<- toGRanges(EnsDb.Hsapiens.v75, feature="gene")

annoData_Mmusculus <<- toGRanges(EnsDb.Mmusculus.v79, feature="gene")

shinyApp(ui = shinyUI(deUI), server = shinyServer(deServer))

