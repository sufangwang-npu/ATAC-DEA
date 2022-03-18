# Build up the ui
# The source file used in this section
options(shiny.maxRequestSize = 30*1024^2)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)

deUI <- function(){

  # Build a drop down menu to notice users
  # Build the head
  header <- dashboardHeader(title = "ATAC-DEA")

  # Build the sidebar
  sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Instruction", tabName="instruction", selected=TRUE,
               icon=icon("glyphicon glyphicon-book", lib="glyphicon")),
      menuItem("Upload Your Data",tabName="dataUpload",
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
                          tabPanel(title = "Introduction",
                                   h1("Introduction"),
                                   p("ATAC-DEA is a web-based platform to help you analyze ATAC-seq data and plot high quality figures.
                                   The ATAC-DEA allows users to visualize differential expression(DE) peaks and do annotation.
                                   Explore the app's features with the example data set pre-loaded.
                                   Upload your genes Expression data first,then submit your data."),
                                   HTML("<p>ATAC-DEA website:
                                        <a href='http://www.atac-dea.xyz:3838/ATAC-DEA'>http://www.atac-dea.xyz:3838/ATAC-DEA</a></p>
                                        <p>ATAC-DEA source:
                                        <a href='https://github.com/sufangwang-npu/ATAC-DEA'>https://github.com/sufangwang-npu/ATAC-DEA</a></p>
                                        <p>Operating system(s): Platform independent</p>
                                        <p>Programming language: R<p>"),
                                   h2("Required files in pretreatment"),
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



                                   h2("Flow Chart"),
                                   p("The flow chart is shown as below"),
                                   tags$img(src = "Flowchart.jpg", width = "747px", height= "528px")

                          ),
                          tabPanel(title = "Data Pretreatment",
                                   h1("ATAC-DEA Data Pretreatment Guide"),
                                   h2("1.	Organize your data files"),
                                   HTML("<p>To do data pretreatment, three files are needed. We recommend organizing all of
                                   the data files in one directory, and placing <b>mapping results files</b> (generally .bam
                                   file) in ‘./reads’; placing <b>peak calling results files</b> (generally .bed file) in ‘./peaks’.
                                        We use some example files in ‘ATAC-DEA/ATAC-DEA_DataPretreatment/extra’ to hint.</p>"),
                                   h2("2. Install dependency packages"),
                                   p("Before the work, some dependency packages needs to be installed"),
                                   HTML('<pre><code class="language-bash">install.packages(c("shiny", "shinyFiles", "stringr"))
                      \nif (!require("BiocManager", quietly = TRUE))
                      \ninstall.packages("BiocManager")
                      \nBiocManager::install("DiffBind")</code></pre>'),
                                   h2("3.	Run pretreatment shiny app"),
                                   HTML('<p>We design a local shiny app to help you do pretreatment. Firstly, you can download the source of ATAC-DEA including
                 Data Pretreatment APP in Github <a href="https://github.com/sufangwang-npu/ATAC-DEA">ATAC-DEA</a>. After that, you can run
                      these commands in terminal to start the app:</p>'),
                                   HTML("<pre><code class='language-bash'>#cd ./ATAC-DEA/ATAC-DEA_DataPretreatment
                      \n#R -e 'shiny::runApp ()'</code></pre>"),
                                   p("After the all starts successfully, type the url (http://127.0.0.1:xxxx typically) in your browser."),
                                   p("Or, you can run this app in Rstudio."),
                                   h2("4.	The panel of ATAC-DEA_DataPretreatment"),
                                   HTML("
                 <p>ATAC-DEA_DataPretreatment helps you build up the sample sheet that DiffBind needs.</p>
                 <p>1)	Open ‘Create DEList’ tab. If you have completed DEList, you can directly upload it and you can also edit it in this panel.</p>
                 <p>2)	Set the work space, if the datapath in your DEList is relative path. The final result file ‘peak_data_collection’ file will also be sived in the work space</p>
                 <p>3)	Input the information of your sample. The DEList in ATAC-DEA/ATAC-DEA_DataPretreatment/extra shows what each variable represents.</p>
                 <p>4)	Select the mapping results file and peak calling results file of each sample. These two buttons only record the relative paths, so don't need to worry about the size of files.</p>
                 <p>5)	Select ‘Insert’ to add a new line and select ‘remove’ to remove the latest line.</p>
                 <p>6)	You can check your DEList here</p>
                 <p>7)	When both DEList and work space are done, you can select ‘complete’ to go to the next step. If one of these two steps has not been done, ‘complete’ button will not response.</p>
                      "),
                                   tags$img(src = "1.jpg", width = "605px", height = "371px"),
                                   p(),
                                   HTML("8) The page will automatically redirect to readCount tab. You can check your DEList again. If there is no mistake, you can press ‘Do Analysis!’ to start counting the reads. It will take some time, depending on your computer."),
                                   tags$img(src = "2.png", width = "610px", height = "279px"),
                                   p(),
                                   HTML("9) When counting is done, a message window will pop out and a file named peak_data_collection will be saved in your work space. Now you can upload this file in ATAC-DEA to explore your data."),
                                   tags$img(src = "3.png", width = "605px", height = "370px"),
                                   ),
                          tabPanel(title = "Tutorial",
                                   h1("ATAC-DEA Quick-start Tutorial"),
                                   p("This quick-start tutorial will guide you through using the example data to do DE & annotation
                                     analysis, and visualize the result."),
                                   p("Open http://www.atac-dea.xyz:3838/ATAC-DEA"),
                                   p("1.	In the ‘Upload Your Data’ tab, you can select to use example or your own data"),
                                   p("1)	Example data: no files need to be uploaded. ATAC-DEA will use data of breast cancer cell
                                     to help you explore it."),
                                   p("2)	Your own data: upload your peak collection result file generated in pretreatment step and
                                     peak calling results files of your data. You can check your data in ‘Check your Sample Sheet’ box."),
                                   HTML("<p>2.	Select the organism of your data: ATAC-DEA now provides ‘<i>Homo sapiens</i>’ and ‘<i>Mus musculus</i>’
                                        to select. When all above is done, select ‘Go Next’ button</p>"),
                                   tags$img(src = "Tutorial1.jpg", width = "633px", height = "390px"),

                                   p("3.	A new sub tab ‘Contrast Design’ will be created and you can design your contrast model.
                                     You can use default option to explore ATAC-DEA. The detail of other options is listed in the paper.
                                     After that, select ‘Do analysis’ to process the analysis and it will take some time."),
                                   p("4.	ATAC-DEA will analysis all possible contrast results in your contrast model and show them at
                                     the bottom of the box for you to select. Because the table is interactive, you can directly select
                                     one for further analysis."),
                                   tags$img(src = "Tutorial2.jpg", width = "633px", height = "390px"),

                                   p("5.	Four new tabs will be shown in sidebar, and you can explore the analysis results by selecting
                                     ‘Data Report’, ‘DE Analysis’ and ‘Peak Annotation’"),
                                   tags$img(src = "Tutorial3.jpg", width = "633px", height = "390px"),

                                   p("6.	Filter is used to set threshold and method (DESeq2 and edgeR are provided) in volcano plot."),
                                   p("7.	Volcano plot is interactive, so you can use brush to select some dots and check their information."),
                                   tags$img(src = "Tutorial4.jpg", width = "633px", height = "390px"),

                                   p("8.	Explore the result of DE analysis in ‘Volcano’, ‘PCA’, ‘Heatmap’ and ‘KEGG&GO’ sub tabs"),
                                   tags$img(src = "Tutorial5.jpg", width = "633px", height = "390px"),

                                   p("9.	In ‘Peak Annotation’ tab, you should select the samples you would like to annotate
                                     (more than 2 and no more than 5 otherwise ATAC-DEA cannot get overlaps)"),
                                   p("10.	ATAC-DEA will show all overlaps according to the samples you selected.
                                     After that select ‘Do Analysis!’ and select the sub tabs to explore the results of annotation."),
                                   tags$img(src = "Tutorial6.jpg", width = "633px", height = "390px"),
                          )
                          )),


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
                                       actionButton("goContrast","Go Next")
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
                               plotOutput("dba_Volcano", brush = "plot_brush")%>% withSpinner(color="#3c8cbc"),
                               downloadButton("downloadVolcanoPlot_DE", "Download Volcano Plot")),

                           box(title="Bar Plot", status = "primary", solidHeader = TRUE,
                               collapsible = TRUE, width = 4,
                               plotOutput("dba_BarPlot")%>% withSpinner(color="#3c8cbc")),

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
                               plotOutput("dba_PCA")%>% withSpinner(color="#3c8cbc"),
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
                               plotOutput("dba_heatmap")%>% withSpinner(color="#3c8cbc"),
                               downloadButton("downloadHeatmapPlot_DE", "Download Heatmap Plot"))
                         )),


                tabPanel("KEGG & GO",
                         fluidRow(
                           tabBox(
                             id = "DEpathway", width = 12,
                             tabPanel("KEGG Enrichment",
                                      plotOutput("DE_KEGGplot")%>% withSpinner(color="#3c8cbc"),
                                      downloadButton("downloadKEGGplot_DE", "Download KEGG")

                                      ),
                             tabPanel("GO Biological Process",
                                      plotOutput("DE_GOBiologicalProcessPathwayPlot")%>% withSpinner(color="#3c8cbc"),
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
