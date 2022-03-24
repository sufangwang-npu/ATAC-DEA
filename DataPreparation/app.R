library(shiny)
library(shinyFiles)
library(stringr)
library(bslib)
library(DiffBind)

ui <- fluidPage(

  theme = bs_theme(version = 4, bootswatch = "cosmo"),

  titlePanel("ATAC-DEA Data Preparation"),

  fluidRow(
      tabsetPanel(id="dataPretreatment",
        tabPanel("Guide", icon = icon("glyphicon glyphicon-book", lib="glyphicon"),
                 fluidRow(
                   column(1),
                   column(10,
                 h1("ATAC-DEA Data Preparation Guide."),
                 h2("1.	Organize your data files"),
                 HTML("<p>To do data preparation, three files are needed. We
                      recommend organizing all of the data files in one directory,
                      and placing <b>mapping results files</b> (generally .bam file) in ‘./reads’;
                      placing <b>peak calling results files</b> (generally .bed file) in ‘./peaks’.
                      We use some example files in ‘ATAC-DEA/DataPreparation/extra’ to hint.</p>"),
                 h2("2. Install dependency packages"),
                 p("Before the work, some dependency packages needs to be installed"),
                 HTML('<pre><code class="language-bash">install.packages(c("shiny", "shinyFiles", "stringr"))
                      \nif (!require("BiocManager", quietly = TRUE))
                      \ninstall.packages("BiocManager")
                      \nBiocManager::install("DiffBind")</code></pre>'),
                 h2("3.	Run preparation shiny app."),
                 HTML('<p>We design a local shiny app to help you do preparation. Firstly, you can download the source of ATAC-DEA including
                 Data Preprartion Application in Github <a href="https://github.com/sufangwang-npu/ATAC-DEA">ATAC-DEA</a>. After that, you can run
                      these commands in terminal to start the app:</p>'),
                 HTML("<pre><code class='language-bash'>#cd ./ATAC-DEA/DataPreparation
                      \n#R -e 'shiny::runApp ()'</code></pre>"),
                 p("After the all starts successfully, type the url (http://127.0.0.1:xxxx) in your browser."),
                 p("Or, you can run this app in Rstudio."),
                 h2("4.	The panel of DataPreparation"),
                 HTML("
                 <p>DataPreparation helps you build up the sample sheet that DiffBind needs.</p>
                 <p>1)	Open ‘Create DEList’ tab. If you have completed DEList, you can directly upload it and you can also edit it in this panel.</p>
                 <p>2)	Set the work space, if the datapath in your DEList is relative path. The final result file ‘peak_data_collection’ file will also be sived in the work space</p>
                 <p>3)	Input the information of your sample. The DEList in ATAC-DEA/DataParation/extra shows what each variable represents.</p>
                 <p>4)	Select the mapping results file and peak calling results file of each sample. These two buttons only record the relative paths, so don't need to worry about the size of files.</p>
                 <p>5)	Select ‘Insert’ to add a new line and select ‘remove’ to remove the latest line.</p>
                 <p>6)	You can check your DEList here</p>
                 <p>7)	When both DEList and work space are done, you can select ‘complete’ to go to the next step. If one of these two steps has not been done, ‘complete’ button will not response.</p>
                      "),
                 tags$img(src = "1.jpg", width = "605px", height = "371px"),
                 p(),
                 HTML("8) The page will automatically redirect to readCount tab. You can check your DEList again. If there is no mistake, you can press ‘Do Analysis!’ to start counting the reads. It will take some time, depending on your computer."),
                 p(),
                 tags$img(src = "2.png", width = "610px", height = "279px"),
                 p(),
                 HTML("9) When counting is done, a message window will pop out and a file named peak_data_collection will be saved in your work space. Now you can upload this file in ATAC-DEA to explore your data."),
                 tags$img(src = "3.png", width = "605px", height = "370px"),
                 )
                 )),
        tabPanel("Create DEList",
                 icon = icon("glyphicon glyphicon-th-list", lib="glyphicon"),
                 fluidRow(
                   column(1),
                   column(3,
                          h3("DEList"),
                          p("    If you have completed DEList or sample sheet, you can directly upload it"),
                          fileInput("userSampleSheet", h4("Choose CSV File"),
                                    multiple = TRUE,
                                    accept = c("text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")),
                          actionButton("submitSampleSheet", "Submit", class = "btn-default",
                                       icon = icon("glyphicon glyphicon-open", lib="glyphicon"),
                                       width = "150px"),
                          hr(),
                          h3("Set Work Space"),
                          p("    We recommend to put '/reads' and '/peaks' in the same directory and set is as work space.
                            peak_data_collection file will also be put in this directory"),
                          p(),
                          shinyDirButton("workSpaceBtn", "Choose Work Space" ,
                                        title = "Please select a dir", buttonType = "default",
                                        class = NULL, style = "width:200px", icon = icon("glyphicon glyphicon-folder-open", lib="glyphicon")),
                          verbatimTextOutput("txt_workSpace")

                          ),
                   column(4,
                          h3("Sample Information"),
                          textInput("tissue", h4("Tissue"),
                                    value = "NA", width = "400px"),
                          textInput("factor", h4("Factor"),
                                    value = "NA", width = "400px"),
                          textInput("condition", h4("Condition"),
                                    value = "NA", width = "400px"),
                          textInput("treatment", h4("Treatment"),
                                    value = "NA", width = "400px"),
                          textInput("replicate", h4("Replicate"),
                                    value = "NA", width = "400px")
                          ),
                   column(4,
                          h3("Pathway of Data Files"),
                          tags$b("bamReads"),
                          p(),
                          shinyFilesButton("bamReadsBtn", "Choose a reads file [.bam]" ,
                                        title = "Please select a file:", multiple = FALSE,
                                        buttonType = "default", class = NULL, icon = icon("glyphicon glyphicon-file", lib="glyphicon"),
                                        style = "width:300px"),
                          verbatimTextOutput("txt_file_bamReads"),
                          p(),
                          tags$b("Peaks"),
                          p(),
                          shinyFilesButton("PeaksBtn", "Choose a peak file" ,
                                        title = "Please select a file:", multiple = FALSE,
                                        buttonType = "default", class = NULL, icon = icon("glyphicon glyphicon-file", lib="glyphicon"),
                                        style = "width:300px"),
                          verbatimTextOutput("txt_file_Peakds"),
                          selectInput("PeakCaller","Peaks Caller",
                                      choices = list(".bed" = "bed", ".narrowPeak" = "narrowPeak",
                                                     ".xls" = "xls"), selected = 1, width = "300px"),
                          hr(),
                          h3("Insert or Remove a Line"),
                          actionButton('insertBtn', HTML("<b>Insert</b>"), class = "btn-primary",
                                       icon = icon("glyphicon glyphicon-plus", lib="glyphicon"),
                                       width = "150px"),
                          p(),
                          actionButton('removeBtn', HTML("<b>Remove</b>"), class = "btn-danger",
                                       icon = icon("glyphicon glyphicon-minus", lib="glyphicon"),
                                       width = "150px"),
                          hr(),
                          actionButton("completeBtn", "Complete!", class = "btn-success",
                                       icon = icon("glyphicon glyphicon-ok", lib="glyphicon"),
                                       width = "150px")

                       )
                 ),
                 column(12, hr()),
                 column(2),
                 column(10,
                        tableOutput("DEList")
                        )

                 ),
        tabPanel("readCount",
                 icon = icon("glyphicon glyphicon-stats", lib="glyphicon"),
                 fluidRow(
                   column(12,
                          tableOutput("dbaSampleSheet"),
                          ),
                   column(12,
                          hr(),
                          actionButton('analysisBtn', "Do Analysis !", class = "btn-primary",
                                       icon = icon("glyphicon glyphicon-play", lib="glyphicon"),
                                       width = "150px"),
                          tags$head(tags$script(src = "message-handler.js")),
                          p(),
                          htmlOutput("finishMessage"))
                 ))
      )
    )
)


server <- function(input, output, session) {

  v <- reactiveValues(DEList = NULL, amReadsDatapath = NULL, PeaksDatapath = NULL,workDatapath = NULL)

  observeEvent(input$submitSampleSheet, {
      req(input$userSampleSheet)
      df <- read.csv(input$userSampleSheet$datapath)
      v$DEList <- df
  })


  observe({
    shinyFileChoose(input, "bamReadsBtn", roots = c(root='.'), session = session)

    if(!is.null(input$bamReadsBtn)){
      file_selected<-parseFilePaths(c(root='.'), input$bamReadsBtn)
      v$bamReadsDatapath <- as.character(file_selected$datapath)
      output$txt_file_bamReads <- renderText(as.character(file_selected$datapath))
    }
  })

  observe({
    shinyFileChoose(input, "PeaksBtn", roots = c(root='.'), session = session)

    if(!is.null(input$PeaksBtn)){
      file_selected <- parseFilePaths(c(root='.'), input$PeaksBtn)
      v$PeaksDatapath <- as.character(file_selected$datapath)
      output$txt_file_Peakds <- renderText(as.character(file_selected$datapath))
    }
  })

  observe({
    shinyDirChoose(input,"workSpaceBtn", roots = c(root='.'), session = session)

    if(!is.null(input$workSpaceBtn)){
      dir_selected <- parseDirPath(c(root='.'), input$workSpaceBtn)
      v$workDatapath <- dir_selected
      output$txt_workSpace <- renderText(dir_selected)
    }

  })


  observeEvent(input$insertBtn,{

    result = tryCatch(
      {
        DEList <- data.frame(v$DEList)
        sampleID <- paste(input$tissue, input$replicate, sep = "_")

        sampleLine <- cbind(sampleID, input$tissue, input$factor, input$condition,
                            input$treatment, input$replicate, v$bamReadsDatapath,
                            v$PeaksDatapath, input$PeakCaller)

        if(nrow(DEList) == 0){
          DEList <- rbind(DEList, sampleLine)
          v$DEList <- DEList

        }else{
          sampleLine <- data.frame(sampleLine)

          lineName <- c("SampleID", "Tissue", "Factor", "Condition","Treatment",
                        "Replicate", "bamReads", "Peaks","PeakCaller")

          names(sampleLine) <- lineName
          names(DEList) <- lineName

          DEList <- rbind(DEList, sampleLine)
          v$DEList <- DEList
        }

      },
      error = function(error_condition) {
        showNotification("Please input datapath_2", type = "error")
        v$DEList <- NULL
      },
      finally = {
        return(NULL)
      })

  })



  observeEvent(input$removeBtn,{
      DEList <- data.frame(v$DEList)
      lastLineNum <- nrow(DEList)
      DEList <- DEList[-lastLineNum,]
      v$DEList <- DEList

  })

  output$DEList <- renderTable({

    if(is.null(v$DEList)){return(NULL)}
    result = tryCatch(
      {
        DEList <- data.frame(v$DEList)
        names(DEList) <- c("SampleID", "Tissue", "Factor", "Condition","Treatment",
                           "Replicate", "bamReads", "Peaks","PeakCaller")

        DEList

      },
      error = function(error_condition) {
        showNotification("Please input datapath_1", type = "error")
        return(NULL)
      })
  })

  dba <- reactiveValues(sampleSheet = NULL)


  observeEvent(input$completeBtn, {

   if(is.null(v$DEList)){
      showNotification("DEList is empty", type = "error")
      return(NULL)
   }else if (unlist(input$workSpaceBtn) == 0){
     showNotification("Work Space is empty", type = "error")
     return(NULL)
   }else{
     dba$sampleSheet <- data.frame(v$DEList)
     updateTabsetPanel(session, "dataPretreatment",selected = "readCount")
  }

  })


  output$dbaSampleSheet <- renderTable(dba$sampleSheet)

  observeEvent(input$analysisBtn, {
    working_directory <- v$workDatapath
    withProgress(message = 'Doing Analysis', value = 0,{
      result = tryCatch(
        {
        setwd(working_directory)
        samples <- data.frame(dba$sampleSheet)
        control <- data.frame("ControlID"=rep(NA, times=nrow(samples)), "bamControl"=rep(NA, times=nrow(samples)))
        df <- cbind(samples[,1:7], control, samples[,8:9])

        incProgress(0.3, detail = "making dba project")
        dba_analysis <- dba(sampleSheet=df)

        incProgress(0.3, detail = "counting")
        dba_analysis <- dba.count(dba_analysis, bUseSummarizeOverlaps = FALSE)
        save(dba_analysis,file=paste("peak_data_collection", sep=""))

        session$sendCustomMessage(type = 'Message',message = 'Successfully！')
        output$finishMessage <- renderText("<b><font color='blue'>Results file has been saved in,now you can go to ATAC-DEA to visualize your data!</font></b>")
        },
        error = function(error_condition) {
          session$sendCustomMessage(type = 'Message',message = 'Oops！something goes wrong!')
          output$finishMessage <- renderText("<b><font color = 'red'>Something is wrong, please check your DEList or data files</font></b>")
        })

    })

  })



}

shinyApp(ui, server)
