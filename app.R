options(shiny.maxRequestSize = 30*1024^2)
# packages used in this app
library(shiny)
library(shinydashboard)
library(shinyjs)
library(BiocManager)
# options(repos = BiocManager::repositories())
library(DiffBind)
library(RColorBrewer)
library(bslib)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(org.Hs.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
# packages used in this app
library(DT)
#source("volcano_plot.R")
source("server.R")
source("ui.R")


annoData_Hsapines <<- toGRanges(EnsDb.Hsapiens.v75, feature="gene")

annoData_Mmusculus <<- toGRanges(EnsDb.Mmusculus.v79, feature="gene")

shinyApp(ui = shinyUI(deUI), server = shinyServer(deServer))

