# Build up the ui
# The source file used in this section
options(shiny.maxRequestSize = 30*1024^2)

# packages used in this app
library(shiny)
library(shinydashboard)
library(BiocManager)
# options(repos = BiocManager::repositories())
library(DiffBind)
library(shinyjs)
library(bslib)
library(ChIPpeakAnno)
library(ggplot2)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(org.Hs.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
# packages used in this app
library(DT)
library(reactome.db)
library(RColorBrewer)

# source("uiBody.R")


deUI <- function(){
  # Build a drop down menu to notice users
  # Build the head
  header <- dashboardHeader(title = "ATAC-DEA")

  # Build the sidebar
  sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Instruction", tabName="instruction",
               icon=icon("glyphicon glyphicon-book", lib="glyphicon")),
      menuItem("Upload Your Data",tabName="dataUpload",selected=TRUE,
               icon=icon("glyphicon glyphicon-open", lib="glyphicon")),
      menuItemOutput("dataReport"),
      menuItemOutput("filter"),
      menuItemOutput("DE_Analysis"),
      menuItemOutput("Peak_Annotation")

    )
  )

  # Build the body
  body <-  dashboardBody(
    tabItems(
      tabItem(tabName = "instruction",
              tabsetPanel(id='instruction',
                          tabPanel(title = "Quick Start",
                                   h1("Get a quick start of ATAC-DEA"),
                                   p("ATAC-DEA is a web-based platform to help you analyze ATAC-seq data and plot high quality figures.
                                   The ATAC-DEA allows users to visualize differential expression(DE) peaks and do annotation.
                                   Explore the app's features with the example data set pre-loaded.
                                   Upload your genes Expression data first,then submit your data."),

                                   h2("Data Requirements"),
                                   p("1. DEList[.csv]:sample sheet contains the main information of each sample"),
                                   p("2. mapping results files[.bam]"),
                                   p("3. peak calling results files[.bed]"),
                                   h3("Example DEList format"),
                                   tags$img(src = "DEList.jpg", width = "633px", height = "160px"),
                                   tags$div(
                                     br(),
                                     HTML("<p><b>SampleID:</b> ID of each sample</p>"),
                                     HTML("<p><b>Tissue/Factor/Condition/Replicate:</b> The variables of each sample</p>"),
                                     HTML("<p><b>bamReads:</b> Path of mapping results files[.bam] (necessary)</p>"),
                                     HTML("<p><b>ControlID:</b> ID of the control data[.bam] (not necessary)</p>"),
                                     HTML("<p><b>bamControl:</b> Path of the control data[.bam] (not necessary)</p>"),
                                     HTML("<p><b>Peaks:</b> Path of peak calling results files[.bed] (necessary)</p>"),
                                     HTML("<p><b>PeakCaller:</b> Suffix of peak calling results files[.bed] (necessary)</p>")),

                                   h2("Data Pretreatment"),
                                   p("The three files need to be processed in order to get the right format which can be operated by
                                     ATAC-DEA. We wrote a script, called readCount.R, for user to easily do pretreatment. It will
                                     take the three needed files and output a peak data collection file which is used as the input
                                     of ATAC-DEA. "),
                                   tags$img(src = "CountR.jpg", width = "613px", height = "340px"),
                                   p("1. change the directory to your own"),
                                   p("2. choose the count option according to your time and computer"),
                                   p("3. save the peak collection result file to your working directory"),

                                   h2("Flow Chart"),
                                   p("The flow chart is shown as below"),
                                   tags$img(src = "Flowchart.jpg", width = "747px", height= "528px")

                          ))),


      tabItem(tabName = "dataUpload",
              # Buile the Page of Data Upload
              tabsetPanel(id="dataUpload_ContrastDesign",
                          tabPanel(title="Data Upload",
                                   fluidRow(
                                     box(
                                       title = "Upload your data", status = "primary", solidHeader = TRUE,
                                       collapsible = TRUE,width=12,
                                       selectInput("data_sel",label="Select the data",choices=c("Example data" ,"Your data")),
                                       fluidRow(
                                         column(6,fileInput("userSample", "Please upload your peak collection result file",multiple = FALSE)),
                                         column(6,fileInput("peakSample", "Please upload your Peak calling results files[.bed]",multiple = TRUE,accept=".bed")),
                                       ),
                                       selectInput("annotationDatabase", label = "Select the organism for annotation",
                                                   choices = c("Homo sapiens", "Mus musculus"),width='50%'),
                                       actionButton("goContrast","Go Next"),
                                     ),
                                     box(
                                       title="Check your Sample Sheet",status="primary",solidHeader=TRUE,
                                       collapsible = TRUE, width=12,
                                       tableOutput("sampleSheet"),
                                     )
                                   ))


              )),

      ######################################
      # Data Report
      tabItem(tabName="dataReport",
              tabsetPanel(
                tabPanel("DE Report",
                         tabPanel("Report",
                                  DT::dataTableOutput("dbaReport"),
                                  downloadButton("downloadReportData", "Download")),
                ),

                tabPanel("DE Annotation",
                         DT::dataTableOutput("dbareport_Annotation_Table"),
                         downloadButton("downloadDE_AnnotationReport", "Download"))

              )

      ),

      #######################################
      tabItem(tabName="DE_Analysis",
              tabsetPanel(
                tabPanel("Volcano",
                         fluidRow(
                           box(title="Description of Volcano Plot",status = "primary",solidHeader = TRUE,
                               collapsible = TRUE,width=12,
                               p("Similar to MA plots, Volcano plots also highlight
                             significantly differentially bound sites and show their
                             fold changes. Here, however, the confidence statistic
                             (FDR or p-value) is shown on a negative log scale,
                             helping visualize the relationship between the magnitude
                             of fold changes and the confidence that sites are
                               differentially bound.")),

                           box(title="Volcano Plot",status = "primary", solidHeader = TRUE,
                               collapsible = TRUE,width=8,
                               plotOutput("dba_Volcano", brush = "plot_brush"),
                               downloadButton("downloadVolcanoPlot_DE", "Download Volcano Plot")),

                           box(title="Bar Plot", status = "primary", solidHeader = TRUE,
                               collapsible = TRUE, width = 4,
                               plotOutput("dba_BarPlot")),

                           box(title="Information",status="primary",solidHeader=TRUE,
                               collapsible= TRUE,width=8,
                               tableOutput("info"))



                         )),


                tabPanel("PCA",
                         fluidRow(
                           box(title="Description of PCA Plot",status = "primary", solidHeader = TRUE,
                               collapsible = TRUE, width=12,
                               p("While the correlation heatmaps already seen are
                             good for showing clustering, plots based on principal
                             components analysis can be used to give a deeper insight
                             into how samples are associated.")),

                           box(title="PCA Plot",status = "primary", solidHeader = TRUE,
                               collapsible = TRUE,width=8,
                               plotOutput("dba_PCA"),
                               downloadButton("downloadPCAPlot_DE", "Download PCA Plot"))
                         )),



                tabPanel("Heatmap",
                         fluidRow(
                           box(title="Description of Heatmap",status = "primary",solidHeader = TRUE,
                               collapsible = TRUE,width=12,
                               p("Another way to view the patterns of binding affinity
                             directly in the differentially bound sites is via a
                             binding affinity heatmap, showing the read scores for
                             some or all of the binding sites.")),

                           box(title="Heatmap",status = "primary", solidHeader = TRUE,
                               collapsible = TRUE,width=8,
                               plotOutput("dba_heatmap"),
                               downloadButton("downloadHeatmapPlot_DE", "Download Heatmap Plot"))
                         )),


                tabPanel("KEGG & GO",
                         fluidRow(
                           tabBox(
                             id = "DEpathway", width = 12,
                             tabPanel("KEGG Enrichment",
                                      plotOutput("DE_KEGGplot"),
                                      downloadButton("downloadKEGGplot_DE", "Download KEGG")

                                      ),
                             tabPanel("GO Biological Process",
                                      plotOutput("DE_GOBiologicalProcessPathwayPlot"),
                                      downloadButton("downloadGOplot_DE", "Download GO")
                                      )
                           ),
                         ))


              )),

      #######################################
      tabItem(tabName="TF_Enrichment",
              tabsetPanel(
                id = "TF_Enrichment",
                tabPanel("Overlap",
                         fluidRow(
                           box(title="Select the Overlap Analysis Object",status="primary",solidHeader = TRUE,
                               collapsible = TRUE,width=8,
                               h4("You are supposed to choose no more than five samples"),
                               uiOutput("overlapSelection"),
                               htmlOutput("overlapSelection_Message"),
                               uiOutput("overlapOption_StartPanel"),
                               uiOutput("startTF_EnrichmentAnalysis_Button"))
                         )
                         ),

                tabPanel(title="Venn Plot", id = "vennPlot",
                         uiOutput("tabPanel_VennPlot")
                         ),

                tabPanel("Genomic Element Distribution",
                         uiOutput("tabPanel_GED")
                         ),

                tabPanel("Overlap Annotate Peaks",
                         uiOutput("tabPanel_AnnotatePeaks")
                         ),

                tabPanel("Individual Annotate Peaks",
                         htmlOutput("tabPanel_IndividualAnnotatePeaks")
                         ),

                tabPanel(title="Overlap Enriched KEGG",
                         uiOutput("tabPanel_KEGGandGO")
                         ),

                tabPanel(title = "Individual Enriched KEGG",
                         htmlOutput("tabPanel_IndividualEnrichedKEGGandGO"))


              )
      )
    )
  )

  # Build ui
  ui <- dashboardPage(
    header,
    sidebar,
    body
  )

}
