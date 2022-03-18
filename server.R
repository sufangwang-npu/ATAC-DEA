
deServer <- function(input, output,session) {

  library(BiocManager)
  library(DiffBind)
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
  library(DT)
  library(reactome.db)
  library(RColorBrewer)
  library(GO.db)


  # DiffBind Reactive Section --------------------------------------------------------
  #
  # Input the dba object after counting
  #
  #####################################
  sampleInput <- reactive({

    if(input$data_sel == "Example data"){
      sample <- load("Data/dba_object")
      tamoxifen_sample <- get(sample)

      return(tamoxifen_sample)

    }else if(input$data_sel == "Your data"){

      req(input$userSample)
      userSample <- load(input$userSample$datapath)
      userSample <- get(userSample)

      return(userSample)

    }
  })
  #################
  #
  # Do Diffbind analysis
  # Return 'dbaobject' which is the core analysis data in DE
  #
  #####################################
  dbaAnalysis <- reactive({
    samples <- sampleInput()
    dbaObject <- dba(samples)

    # Normalizing the data
    dbaObject <- dba.normalize(dbaObject, method = DBA_DESEQ2,
                               normalize = DBA_NORM_DEFAULT, library = DBA_LIBSIZE_DEFAULT,
                               background = FALSE, bRetrieve=FALSE)

    # The first Mode: Specific contrast, disign and contrast is present
    if(input$contrastMode == "Set up a specific contrast"){
      if(input$factor == "Single Factor"){

        # Establishing a model design and contrast
        # Performing the differential analysis
        contrastDesign <- contrastDesignWithThreeLength() # Return the three length vector used as contrast parameter in dba.contrast
        dbaObject <- dba.contrast(dbaObject,contrast=contrastDesign)
        dbaObject <- dba.analyze(dbaObject,bBlacklist=FALSE,bGreylist=FALSE)

      }else if(input$factor == "Multiple Factor"){

        multipleFactorDesign <- vector()
        for(i in input$multipleFactorSelectBox){multipleFactorDesign <- paste(multipleFactorDesign,i,sep=" + ")}
        multipleFactorDesign <- paste("~",multipleFactorDesign,sep="")
        dbaObject <- dba.contrast(dbaObject,design=multipleFactorDesign)
        dbaObject <- dba.analyze(dbaObject,method=DBA_ALL_METHODS,bBlacklist=FALSE,bGreylist=FALSE)

      }}

    # The second Mode: Automatically generate contrast
    if(input$contrastMode == "Set up a all possible contrasts"){

      dbaObject <- dba.contrast(dbaObject)
      dbaObject <- dba.analyze(dbaObject,method=DBA_ALL_METHODS,bBlacklist=FALSE,bGreylist=FALSE)

    }

    return(dbaObject)

  })

  # Build Up contrast like c('Tissue','Responsive','Resistant')
  contrastDesignWithThreeLength <- reactive({

    contrastDesign <- c()
    contrastDesign <- c(input$singleFactorConstractModel,
                        input$contrast_value1,
                        input$contrast_value2)

    return(contrastDesign)

  })

  #  As user press 'DO Analysis' the analysis process will begin
  observeEvent(input$start,{

    withProgress(message = 'Doing Analysis', value = 0.5,{dbaAnalysis()})

  })

  #################
  #
  # contrastSelect()
  #
  #####################################
  contrastSelect <- reactive({

    contrastSelect <- input$allPossibleContrastInModeII_rows_selected

    if(is.null(contrastSelect)){

      contrastSelect <- 1

    }else{

      contrastSelect <- input$allPossibleContrastInModeII_rows_selected

    }

    return(contrastSelect)

  })

  #################


  # ChIPpeakAnno Reactive Section --------------------------------------------------------
  #
  # Integrate ChIPpeakAnno data
  # peakData() is used to read data and return a summary list including 'peakfileannotion' list
  #
  #####################################
  peakData <- reactive({

    if(input$data_sel == "Example data"){

      peakdata <- read.csv("Data/peak_path.csv")
      peakDataFile <- list()
      ChIPpeakAnno_Grange <- list()
      ChIPpeakAnno_Replicate <-list()
      peakName <- vector()

      for(i in 1:length(peakdata$datapath)){

        name <- paste("peak",i,sep="")
        peakDataFile <- c(peakDataFile,name)
        ChIPpeakAnno_Grange <- c(ChIPpeakAnno_Grange,name)
        peakName <- append(peakName,peakDataFile[[i]])
        peakDataFile[[i]] <- file(peakdata$datapath[i])
        ChIPpeakAnno_Grange[[i]] <- toGRanges(peakDataFile[[i]], format="BED", header=FALSE)
        ChIPpeakAnno_Replicate <- c(ChIPpeakAnno_Replicate,ChIPpeakAnno_Grange[[i]])

      }

      names(ChIPpeakAnno_Grange) <- peakdata$name
      names(ChIPpeakAnno_Replicate) <- peakdata$name

      exampleSummary <<- list(peakDataFile = peakDataFile,
                              ChIPpeakAnno_Grange = ChIPpeakAnno_Grange,
                              ChIPpeakAnno_Replicate = ChIPpeakAnno_Replicate,
                              peakName = peakdata$name)

      return(exampleSummary)

    }else if(input$data_sel == "Your data"){

      # User upload the data
      peakDataFile <- list()
      ChIPpeakAnno_Grange <- list()
      ChIPpeakAnno_Replicate <-list()
      peakName <- vector()


      for(i in 1:length(input$peakSample$datapath)){

        name <- paste("peak",i,sep="")
        peakDataFile <- c(peakDataFile,name)
        ChIPpeakAnno_Grange <- c(ChIPpeakAnno_Grange,name)
        peakName <- append(peakName,peakDataFile[[i]])
        peakDataFile[[i]] <- file(input$peakSample$datapath[i])
        ChIPpeakAnno_Grange[[i]] <- toGRanges(peakDataFile[[i]], format="BED", header=FALSE)
        ChIPpeakAnno_Replicate <- c(ChIPpeakAnno_Replicate,ChIPpeakAnno_Grange[[i]])

      }

      names(ChIPpeakAnno_Grange) <- input$peakSample$name
      names(ChIPpeakAnno_Replicate) <- input$peakSample$name

      userSummary <<- list(peakDataFile = peakDataFile,
                           ChIPpeakAnno_Grange = ChIPpeakAnno_Grange,
                           ChIPpeakAnno_Replicate = ChIPpeakAnno_Replicate,
                           peakName = input$peakSample$name)

      return(userSummary)

    }
  })

  #################
  #
  # Get Annotation Databases prepared
  # Annotation database variables are used in global environment
  #
  #####################################

  observeEvent(input$goContrast,{

    withProgress(message = 'Doing Analysis', value = 0.5,{
      if(input$annotationDatabase == 'Homo sapiens'){

        annoData_Hsapines <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
        # annoData_Hsapines <- load("Data/AnnoDatabase/annoData_Hsapines")
        # annoData_Hsapines <- get(annoData_Hsapines)

        annoData <<- annoData_Hsapines
        annoEnsDb <<- EnsDb.Hsapiens.v75
        annoTxDb <<- TxDb.Hsapiens.UCSC.hg19.knownGene
        annoOrg <<- "org.Hs.eg.db"

      }else if (input$annotationDatabase == 'Mus musculus'){

        # annoData_Mmusculus <- toGRanges(EnsDb.Mmusculus.v79, feature="gene")
        annoData_Mmusculus <- load("Data/AnnoDatabase/annoData_Mmusculus")
        annoData_Mmusculus <- get(annoData_Mmusculus)

        annoData <<- annoData_Mmusculus
        annoEnsDb <<- EnsDb.Mmusculus.v79
        annoTxDb <<- TxDb.Mmusculus.UCSC.mm10.knownGene
        annoOrg <<- "org.Mm.eg.db"

      }
    })
  })

  #################
  #
  # Commpute the overlaps of selected peaks
  #
  #####################################
  # Return the overlap
  generateOverlapObject <- reactive({

    peakSummary <- peakData()
    overlaps <- findOverlapsOfPeaks(peakSummary$ChIPpeakAnno_Grange[input$overlapSampleSelection])
    overlaps <-addMetadata(overlaps, colNames="score", FUN=mean)

    return(overlaps)
  })

  #################
  #
  # Peak Annotation
  # ChIPpeakAnno_PeakAnnotation()
  #
  #####################################
  ChIPpeakAnno_PeakAnnotation <- reactive({

    overlaps <- generateOverlapObject()
    overlaps_plot <- overlaps$peaklist[[input$overlapSelection]]

    overlaps.anno <- annotatePeakInBatch(overlaps_plot,AnnotationData=annoData)

    overlaps.anno <- addGeneIDs(overlaps.anno, annoOrg)
    return(overlaps.anno)

  })
  #################
  #
  # KEGG and GO Annotation: output the pathway data table
  #
  #####################################
  DE_Annotation_KEGG <- reactive({

    DEAnnotation_Report_2 <- DE_Annotation()

    DE_KEGG <- getEnrichedPATH(DEAnnotation_Report_2, annoOrg, "reactome.db", maxP= 0.05)

    return(DE_KEGG)

  })


  DE_Annotation_GO <- reactive({

    DEAnnotation_Report_2 <- DE_Annotation()

    DE_GO <- getEnrichedGO(DEAnnotation_Report_2, orgAnn=annoOrg, condense=TRUE,maxP = 0.05)

    return(DE_GO)

  })

  # Function to produce the KEGG pathway
  ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table <- reactive({

    overlaps.anno <- ChIPpeakAnno_PeakAnnotation()

    withProgress(message = 'Making plot', value = 0.3,{

      BiologicalPath <- getEnrichedPATH(overlaps.anno, annoOrg, "reactome.db", maxP = 0.05)

      return(BiologicalPath)

    })
  })

  ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table <- reactive({

    overlaps.anno <- ChIPpeakAnno_PeakAnnotation()

    withProgress(message = 'Making plot', value = 0.3,{

      BiologicalPath <- getEnrichedGO(overlaps.anno, orgAnn=annoOrg, condense=TRUE,maxP = 0.05)

      return(BiologicalPath)

    })
  })


  # Function to produce the KEGG pathway
  ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i <- function(SampleNum){

    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput <- overlapSelection_UserInput[SampleNum]
    peakSummary <- peakData()

    Sample.anno <- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[SampleNum]],AnnotationData=annoData)

    BiologicalPath <- getEnrichedPATH(Sample.anno, annoOrg, "reactome.db", maxP= 0.05)

    return(BiologicalPath)

  }


  # Function to produce the GO bp pathway
  ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i <- function(SampleNum){

    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput <- overlapSelection_UserInput[SampleNum]
    peakSummary <- peakData()

    Sample.anno <- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[SampleNum]],AnnotationData=annoData)

    BiologicalPath <- getEnrichedGO(Sample.anno, orgAnn=annoOrg, condense=TRUE,maxP = 0.05)

    return(BiologicalPath)

  }
  #################



  # UI Section --------------------------------------------------------------
  #
  # Data Upload
  #
  #####################################

  output$userSample <- renderUI({

    if(input$data_sel == "Your data"){

      #Upload user's data
      fileInput("userSample", "Please upload your dba object after contrast",multiple = FALSE)

    }
  })

  #################
  #
  # ATAC-DEA Menu
  #
  #####################################
  # Side bar menu output: Data Report
  output$dataReport <- renderMenu({

    validate(
      need(input$start != "", "")
    )
    if(input$start!=0){
      menuItem("Data Report",tabName="dataReport",icon=icon("glyphicon glyphicon-list-alt",lib="glyphicon"))
    }

  })

  # Side bar menu output: DE analysis
  output$DE_Analysis <- renderMenu({

    validate(
      need(input$start != "", "")
    )
    if(input$start!=0){

      menuItem("DE Analysis",tabName="DE_Analysis",icon=icon("glyphicon glyphicon-equalizer",lib="glyphicon"))

    }

  })

  # Side bar menu output: TF Enrichment
  output$Peak_Annotation <- renderMenu({

    validate(
      need(input$start != "", "")
    )
    if(input$start!=0){

      menuItem("Peak Annotation",tabName="TF_Enrichment",icon=icon("glyphicon glyphicon-stats",lib="glyphicon"))

    }

  })

  # Side bar menu output: Control
  output$filter <- renderMenu({

    validate(
      need(input$start != "", "")
    )
    if(input$start!=0){

      menuItem("Filter",tabName="Filter",icon=icon("glyphicon glyphicon-filter",lib="glyphicon"),
               sliderInput("fdr",label="FDR",min=1e-10, max=0.05, value=0.005),
               sliderInput("fold_change",label="log2FoldChange",min=0.5,max=5,value=0.5),
               selectInput("plot_method",label="Analysis Method",choices=c("DESeq2"=DBA_DESEQ2,"EdgeR"=DBA_EDGER)))

    }
  })

  #################
  #
  # Contrast Model Selection Panel UI
  #
  #####################################

  # Contrast panel Design
  # When user successfully upload the data and press 'Go Next' button,
  # The contrast design Panel will show
  observeEvent(input$goContrast,{
      insertTab(
        inputId = "dataUpload_ContrastDesign",
        tabPanel(title="Contrast Design", id="dbaContrastDesign",
                 fluidRow(
                   box(
                     title="Design your contrast", status = "primary", solidHeader = TRUE,
                     collapsible = TRUE, width = 12,
                     column(12,selectInput("contrastMode","Which kind of contrast mode would you like to use?",
                                           choices = c("Set up a specific contrast",
                                                       "Set up a all possible contrasts"),
                                           selected = 1,width = '80%')),
                     column(12,uiOutput("numberOfFactor")),
                     column(12,uiOutput("constractDesignSelection")),
                     column(8,actionButton("start","Do Contrast")),
                     column(12,uiOutput("contrastResultReturn")),
                     column(12,DT::dataTableOutput("allPossibleContrastInModeII"))

                   ),
                 )),
        target = "Data Upload",
        position = "after"
      )

    updateTabsetPanel(
      session = getDefaultReactiveDomain(),
      inputId = "dataUpload_ContrastDesign",
      selected = "Contrast Design"

    )


  })

  # Dynamic ui:contrast model
  # value 1 in contrast
  output$contrastModel_value1 <- renderUI({

    sampleSheet <- sampleInput()
    sampleSheet <- sampleSheet$samples
    singleFactorConstractModel_UserInput <- input$singleFactorConstractModel
    condition <- unique(sampleSheet[singleFactorConstractModel_UserInput])
    condition <- as.vector(condition[[1]])

    selectInput("contrast_value1",label="Contrast Model",
                choices=condition,width = '80%')

  })

  # value 2 in contrast
  output$contrastModel_value2 <- renderUI({

    sampleSheet <- sampleInput()
    sampleSheet <- sampleSheet$samples
    singleFactorConstractModel_UserInput <- input$singleFactorConstractModel
    condition <- unique(sampleSheet[singleFactorConstractModel_UserInput])
    condition <- as.vector(condition[[1]])
    deleteRepetitiveValue <- which(condition==input$contrast_value1) # Delete value in contrast 1
    conditionAfterDeletion <- condition[-deleteRepetitiveValue]

    selectInput("contrast_value2",label="Contrast Model",
                choices=conditionAfterDeletion,width = '80%')

  })

  # Return corresponding contrast ui in Mode I
  output$constractDesignSelection <- renderUI({
    if(input$contrastMode == "Set up a specific contrast"){
      uiOutput("multipleFactorPanel")
    }else if(input$contrastMode == "Set up a all possible contrasts"){
      return()
    }
  })

  output$numberOfFactor <- renderUI({
    if(input$contrastMode == "Set up a specific contrast"){
      selectInput("factor","Single or Multiple Factor analysis",
                  choices=c("Single Factor","Multiple Factor"),selected = 1,width = '80%')
    }
  })


  # Return corresponding contrast ui in ModeII
  observeEvent(input$start,{

    # The data table of Contrast Result
    output$allPossibleContrastInModeII <- DT::renderDataTable({
      dbaObject <- dbaAnalysis()
      dba.show(dbaObject,bContrasts=TRUE)
    }, options = list(dom = 't',ordering=FALSE,autoWidth=TRUE),
    selection=list(mode="single"),

    )

    output$contrastResultReturn <- renderUI({

      column(12,h4("Contrast Result"),p("You choose the contrast to be analyzed"))

    })

    # When user press 'Do Contrast', the button will be replaces to 'Do Analysis'
    # whose label is input$start
    updateActionButton(
      session,
      "start",
      label = "Do Analysis",
    )

    # When user press 'Do Analysis, the button will be disable
  })

  #################
  #
  # Multi-factor Contrast Model Establishment Panel UI
  #
  #####################################
  output$multipleFactorPanel <- renderUI({

    sampleSheet <- sampleInput()
    sampleSheet <- sampleSheet$samples
    # Create a contrast object data.frame
    sampleSheet_Tissues <- sampleSheet$Tissue
    sampleSheet_Factor <- sampleSheet$Factor
    sampleSheet_Condition <- sampleSheet$Condition
    sampleSheet_Treatment <- sampleSheet$Treatment
    sampleSheet_Replicate <- sampleSheet$Replicate
    sampleSheet_Caller <- sampleSheet$PeakCaller

    sampleSheetContrastObject <-data.frame(Tissue=sampleSheet_Tissues,Factor=sampleSheet_Factor,
                                           Condition=sampleSheet_Condition,Treatment=sampleSheet_Treatment,
                                           Replicate=sampleSheet_Replicate,Caller=sampleSheet_Caller)

    sampleSheetContrastObject_Point <- 1 # A point to find the name of corresponding column

    contrastFactorSelection <- vector()

    for(i in sampleSheetContrastObject){

      temp <- i[1]

      if(!all(i == temp)){

        contrastFactorSelection <- cbind(contrastFactorSelection,
                                         names(sampleSheetContrastObject[sampleSheetContrastObject_Point]))
      }

      sampleSheetContrastObject_Point <- sampleSheetContrastObject_Point + 1

    }



    if(input$factor == "Multiple Factor"){

      # Multiple Factor
      box(title="Select the Multiple Factors",
          width=8,
          checkboxGroupInput("multipleFactorSelectBox","Choose the factors",
                             choices = contrastFactorSelection,
                             selected = 2)
      )

    }else if(input$factor == "Single Factor"){

      # Single Factor panel
      box(title="Select the Factors",width=8,
          selectInput("singleFactorConstractModel","Choose the contrast factor",
                      choices = contrastFactorSelection,selected=1,width = '80%'),
          uiOutput("contrastModel_value1"),
          uiOutput("contrastModel_value2")
      )

    }
  })

  #################
  #
  # Dynamic peak selection panel UI
  #
  #####################################

  # Overlap Selection
  output$overlapSelection <- renderUI({

    peakSummary <- peakData()
    peakName <- peakSummary$peakName

    checkboxGroupInput("overlapSampleSelection",label="Select the sample",
                       choices = peakName,
                       selected = 2)

  })

  # return the message according to user's input
  overlapSelection_Return <- reactive({

    overlapSelection_UserInput <- input$overlapSampleSelection

    if(length(overlapSelection_UserInput)<2){

      message <- "Please choose more samples."
      userInput <- NULL
      resultReturn <- list(message=message,userInput=userInput)

      return(resultReturn)

    }else if(length(overlapSelection_UserInput)>=2 & length(overlapSelection_UserInput)<=5){

      message <- "The number of samples is alright."
      userInput <- overlapSelection_UserInput
      resultReturn <- list(message=message,userInput=userInput)

      return(resultReturn)

    }else if(length(overlapSelection_UserInput)>5){

      message <- "There are too many samples"
      userInput <- NULL
      resultReturn <- list(message=message,userInput=userInput)

      return(resultReturn)

    }
  })

  output$overlapSelection_Message <- renderText({

    overlapSelection_UserInput <- overlapSelection_Return()
    overlapSelection_UserInput_Message <- overlapSelection_UserInput[["message"]]

    if(is.null(overlapSelection_UserInput[["userInput"]])){
      paste("<p style='color:red'><b>",overlapSelection_UserInput_Message,"</b></p>")

    }else{
      paste("<p style='color:blue'><b>",overlapSelection_UserInput_Message,"</b></p>")
    }

  })

  # When the sample is no more than 5, the button will appear
  output$startTF_EnrichmentAnalysis_Button <- renderUI({

    overlapSelection_UserInput <- input$overlapSampleSelection

    if(length(overlapSelection_UserInput)>=2 & length(overlapSelection_UserInput)<=5){

      actionButton("startTF_EnrichmentAnalysis",tags$b("Do Analysis!"),class="btn-primary",style="color:white")

    }else{ return() }

  })

  output$overlapOption_StartPanel <- renderUI({

    overlapSelection_UserInput <- input$overlapSampleSelection

    if(length(overlapSelection_UserInput)>=2 & length(overlapSelection_UserInput)<=5){

      peakNameOverlap <- createOverlapOption()
      selectInput("overlapSelection", h4("Select overlap sample"),
                  choices = peakNameOverlap, selected = 1 )

    }else{ return() }
  })

  #################
  #
  # Overlap Option UI: This option will return the combination of different overlaps
  #
  #####################################
  createOverlapOption <- reactive({

    temp <- vector()
    overlaps <- generateOverlapObject()
    peakNameOverlap <- names(overlaps$peaklist)

    for(i in 1:length(peakNameOverlap)){

      if(grepl("///",peakNameOverlap[i])){
        temp <- c(temp,peakNameOverlap[i])

      }
    }

    peakNameOverlap<-vector()
    peakNameOverlap <- temp
    return(peakNameOverlap)

  })

  output$overlapOption <- renderUI({

    peakNameOverlap <- createOverlapOption()
    selectInput("overlapSelection_Plot", h4("Select overlap sample"),
                choices = peakNameOverlap, selected = 1 )

  })

  ################
  #
  # GED plot Panel UI
  #
  #####################################
  output$tabPanel_GED <- renderUI({

    input$startTF_EnrichmentAnalysis

    if(is.null(input$startTF_EnrichmentAnalysis)){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please select the sample in 'Overlap' Panel !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis == 0){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please press the 'Do Analysis' button !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis >= 1){

      fluidRow(

        box(title="Genomic Element Distribution of Duplicates", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width=12,
            plotOutput("ChIPpeakAnno_DuplicateDistribution")%>% withSpinner(color="#3c8cbc"),
            actionButton("duplicateDistributionPlot","Plot",class="btn-primary",style="color:white"),
            downloadButton("downloadDuplicateDistribution_AN", "Download GEDD Plot")),

        box(title="Genomic Element Distribution of Overlaps", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width=6,
            plotOutput("ChIPpeakAnno_OverlapDistribution")%>% withSpinner(color="#3c8cbc"),
            actionButton("overlapDistributionPlot","Plot",class="btn-primary",style="color:white"),
            downloadButton("downloadOverlapDistributionPlot_AN", "Download GEDO Plot"))
      )
    }
  })

  #################
  #
  # Venn Plot Panel UI
  #
  #####################################
  output$tabPanel_VennPlot <- renderUI({

    input$startTF_EnrichmentAnalysis

    isolate({
      if(is.null(input$startTF_EnrichmentAnalysis)){
        fluidRow(
          box(title="Warring Message", status = "warning", solidHeader = TRUE,
              collapsible = TRUE, width=8,
              br(),
              tags$p(style='color:red',"Please select the sample in 'Overlap' Panel !"),
              br(),

          ))
      }else if(input$startTF_EnrichmentAnalysis == 0){
        fluidRow(
          box(title="Warring Message", status = "warning", solidHeader = TRUE,
              collapsible = TRUE, width=8,
              br(),
              tags$p(style='color:red',"Please press the 'Do Analysis' button !"),
              br(),

          ))
      }else if(input$startTF_EnrichmentAnalysis >= 1){
        fluidRow(
          box(title="Venn Plot", status = "primary", solidHeader = TRUE,
              collapsible = TRUE,width=8,
              plotOutput("ChIPpeakAnno_Venn")%>% withSpinner(color="#3c8cbc"),
              downloadButton("downloadVennPlot_AN", "Download Venn Plot")),

          tabBox(title ="Option", width=4,id = "overlapOption", height = "250px",
                 tabPanel("Venn",
                          uiOutput("vennPlotOption_SelectColor"),
                          actionButton("vennPlot",tags$b("Plot"),class="btn-primary",style="color:white"))

          ))
      }
    })
  })

  #################
  #
  # Venn Plot Color Selection Box UI
  #
  #####################################
  output$vennPlotOption_SelectColor <- renderUI({
    selectInput("selectColor", h4("Select Color"),
                choices = list("Color 1" = "color1",
                               "Color 2" = "color2",
                               "Color 3" = "color3"),
                selected = 1)
  })

  vennPlotOption_Color <- reactive({
    vennPlotOption_UserInput <- input$overlapSampleSelection
    vennPlotOption_UserInput_Length <- length(vennPlotOption_UserInput)
    if(vennPlotOption_UserInput_Length == 2){

      switch(input$selectColor,
             color1 = c("#8dd3c7","#bebada"),
             color2 = c("#66c2a5","#fc8d62"),
             color3 = c("#a6cee3","#b2df8a"))
    }
    else if(vennPlotOption_UserInput_Length == 3){

      switch(input$selectColor,
             color1 = c("#8dd3c7","#fdb462","#bebada"),
             color2 = c("#66c2a5","#fc8d62","#8da0cb"),
             color3 = c("#a6cee3","#1f78b4","#b2df8a"))
    }
    else if(vennPlotOption_UserInput_Length == 4){

      switch(input$selectColor,
             color1 = c("#8dd3c7","#fdb462","#bebada","#fb8072"),
             color2 = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"),
             color3 = c("#a6cee3","#1f78b4","#b2df8a","#33a02c"))
    }
    else if(vennPlotOption_UserInput_Length == 5){

      switch(input$selectColor,
             color1 = c("#8dd3c7","#fdb462","#bebada","#fb8072","#80b1d3"),
             color2 = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854"),
             color3 = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99"))
    }

  })

  #################
  #
  # Overlap Peak Annotation: Table and Pie Plot
  #
  #####################################
  output$tabPanel_AnnotatePeaks <- renderUI({

    input$startTF_EnrichmentAnalysis

    if(is.null(input$startTF_EnrichmentAnalysis)){
      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please select the sample in 'Overlap' Panel !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis == 0){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please press the 'Do Analysis' button !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis >= 1){

      fluidRow(
        box(title="Annotate the Overpal Peaks to the Promoter Regions Reporter", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            verbatimTextOutput("ChIPpeakAnno_PeakAnnotation_Table")%>% withSpinner(color="#3c8cbc"),
            downloadButton("downloadAnnotatePeaksData", "Download Annotation Table")),

        box(title="Overlap Pie Plot",status = "primary", solidHeader = TRUE,
            collapsible = TRUE,width=4,
            plotOutput("ChIPpeakAnno_PeakAnnotation_Pie")%>% withSpinner(color="#3c8cbc"),
            downloadButton("downloadAnnotatePiePlot", "Download the Overlap Pie plot"))
      )
    }
  })

  #################
  #
  # Individual Peak Annotation: Table and Pie Plot
  #
  #####################################
  output$tabPanel_IndividualAnnotatePeaks <- renderPrint({

    overlapSelection_UserInput <- input$overlapSampleSelection
    input$startTF_EnrichmentAnalysis

    if(is.null(input$startTF_EnrichmentAnalysis)){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please select the sample in 'Overlap' Panel !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis == 0){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please press the 'Do Analysis' button !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis >= 1){

      fluidPage(
        for(i in 1:length(overlapSelection_UserInput)){

          print(
            box(title=paste("Annotate the", overlapSelection_UserInput[i] ,"Peaks to the Promoter Regions Reporter",sep=" "),
                status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width=8,
                verbatimTextOutput(paste("ChIPpeakAnno_PeakAnnotation_Table",i,sep = "_")),
                downloadButton(paste("downloadAnnotatePeaksData",i,sep="_"),
                               paste("Download the",overlapSelection_UserInput[i], "Annotation Table"))
            )
          )
          print(box(title=paste(overlapSelection_UserInput[i],"Pie Plot",sep=" "),status = "primary", solidHeader = TRUE,
                    collapsible = TRUE,width=4,
                    plotOutput(paste("ChIPpeakAnno_PeakAnnotation_Pie",i,sep="_")),
                    downloadButton(paste("downloadAnnotatePiePlot",i,sep="_"),
                                   paste("Download the",overlapSelection_UserInput[i], "Pie plot"))
            )
          )
        }
      )
    }
  })

  #################
  #
  # Overlap KEGG & GO Plot Panel UI
  #
  #####################################
  # KEGG&GO plot
  output$tabPanel_KEGGandGO <- renderUI({

    input$startTF_EnrichmentAnalysis

    if(is.null(input$startTF_EnrichmentAnalysis)){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please select the sample in 'Overlap' Panel !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis == 0){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please press the 'Do Analysis' button !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis >= 1){

      fluidRow(
        tabBox(
          title = "Overlap Obtain Enriched KEGG & GO", side = "right",
          id = "DEpathway", width = 12,
          tabPanel("KEGG",
                   plotOutput("ChIPpeakAnno_PeakAnnotation_KEGGPathwayPlot"),
                   downloadButton("downloadKEGGplot", "Download GO Plot")

          ),
          tabPanel("GO Biological Process",
                   plotOutput("ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot"),
                   downloadButton("downloadGOBiologicalProcessplot", "Download GO")
          )
        )
      )
    }
  })

  #################
  #
  # Individual KEGG & GO Plot Panel UI
  #
  #####################################
  # Individual Go plot
  output$tabPanel_IndividualEnrichedKEGGandGO <- renderPrint({

    overlapSelection_UserInput <- input$overlapSampleSelection
    input$startTF_EnrichmentAnalysis

    if(is.null(input$startTF_EnrichmentAnalysis)){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please select the sample in 'Overlap' Panel !"),
            br(),

        ))

    }else if(input$startTF_EnrichmentAnalysis == 0){

      fluidRow(
        box(title="Warring Message", status = "warning", solidHeader = TRUE,
            collapsible = TRUE, width=8,
            br(),
            tags$p(style='color:red',"Please press the 'Do Analysis' button !"),
            br(),

        ))
    }else if(input$startTF_EnrichmentAnalysis >= 1){

      fluidRow(

        for(i in 1:length(overlapSelection_UserInput)){

          print(

            tabBox(
              title = paste("Obtain",overlapSelection_UserInput[i] ,"Enriched KEGG & GO"),
              side = "right", id = "IndividualPathway", width = 12,

              tabPanel("KEGG",
                       plotOutput(paste("ChIPpeakAnno_PeakAnnotation_PathwayPlot",i,sep="_")),
                       downloadButton(paste("downloadKEGGplot",i,sep="_"),
                                      paste("Download",overlapSelection_UserInput[i], "KEGG Plot",sep=" "))

              ),

              tabPanel("GO Biological Process",
                       plotOutput(paste("ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot", i, sep = "_")),
                       downloadButton(paste("downloadGOBiologicalProcessplot", i, sep = "_"),
                                      paste("Download",overlapSelection_UserInput[i], "GO Plot",sep=" "))
              )

            )
            )
        }
      )
    }
  })

  #################


  # Table Section -----------------------------------------------------------
  #
  # SampleSheet Table
  #
  #####################################

  output$sampleSheet <- renderTable({
    sampleSheet_list <- sampleInput()
    sampleSheet <- sampleSheet_list$samples

    sampleSheet <- cbind(sampleSheet[,1:6],sampleSheet[c("ControlID")],sampleSheet[,10:11])
    colnames(sampleSheet)[9]<-"DataFormat"
    sampleSheet

  })

  #################
  #
  # Peaks Information Table of Volcano Plot
  #
  #####################################

  dbareport_Annotation_Table_info <- reactive({
    dbaObject <- dbaAnalysis()
    contrastSelect <- contrastSelect()
    AllPeaks_Report_noAnno_BeforeEdit <- dba.report(dbaObject,contrastSelect,bUsePval=TRUE,
                                                    th=100,bNormalized=TRUE,precision=0,
                                                    DataType=DBA_DATA_FRAME)
    AllPeaks_Report_AfterEdit <- addChr2DEReport(AllPeaks_Report_noAnno_BeforeEdit)
    # DEAnnotation_Report is the dbareport object (GRangs)

    DEAnnotation_Report_Table <- as.data.frame(AllPeaks_Report_AfterEdit)

    DEAnnotation_Report_1 <- annotatePeakInBatch(AllPeaks_Report_AfterEdit, AnnotationData=annoData)
    DEAnnotation_Report_2 <- addGeneIDs(DEAnnotation_Report_1, annoOrg)
    DEAnnotation_Report_3 <- as.data.frame(DEAnnotation_Report_2@elementMetadata)
    DEAnnotation_Report_3_Table <- as.data.frame(DEAnnotation_Report_3)

    if(!is.null(!duplicated(DEAnnotation_Report_3_Table$peak))){

      DEAnnotation_Report_3_Table <- DEAnnotation_Report_3_Table[!duplicated(DEAnnotation_Report_3_Table$peak),]

    }


    DEAnnotation.info <- DEAnnotation_Report_3_Table

    DEAnnotation.info

  })

  output$info <- renderTable({

    dbaObject <- dbaAnalysis()
    volcanoPlotDataList <- volcano_plot_MatrixCreate_Plot()
    DEAnnotation.info <- dbareport_Annotation_Table_info()
    volcano_plot_matrix <- volcanoPlotDataList$volcano_plot_matrix

    volcano_plot_matrix <- cbind(volcano_plot_matrix, Symbol=DEAnnotation.info$symbol)

    brushedPoints(volcano_plot_matrix, input$plot_brush, xvar ="log2FC", yvar="log10FDR")

  })

  #################
  #
  # DE Report Table
  #
  #####################################
  # Retrieving the differentially bound sites
  dbareport <- reactive({

    withProgress(message = 'Generating the report', value = 0.5,{

      dbaObject <- dbaAnalysis()
      contrastSelect <- contrastSelect()
      #str(dbaObject)
      tamoxifen.DB <- dba.report(dbaObject,contrast=contrastSelect, th=input$fdr, bUsePval=FALSE,
                                 fold=input$fold_change,
                                 precision=3:5,
                                 file="tamoxifen_report.csv",DataType=DBA_DATA_FRAME)

    })

    tamoxifen.DB ##data <- round(data,digits=2)

  })

  # Report
  output$dbaReport <- DT::renderDataTable({

    dbareport()

  })

  #################
  #
  # DE Annotation Table
  #
  #####################################

  dbareport_Annotation <- reactive({

    withProgress(message = 'Generating the report', value = 0.5,{

      dbaObject <- dbaAnalysis()
      contrastSelect <- contrastSelect()

      tamoxifen.DB <- dba.report(dbaObject,contrast=contrastSelect, th=input$fdr, bUsePval=FALSE,
                                 fold=input$fold_change,
                                 precision=3:5,
                                 file="tamoxifen_report.csv")

    })
    tamoxifen.DB
  })

  DE_Annotation <- reactive({

    DEAnnotation_Report_noAnno <- dbareport_Annotation()
    DEAnnotation_Report <- addChr2DEReport(DEAnnotation_Report_noAnno)
    # DEAnnotation_Report is the dbareport object (GRangs)

    DEAnnotation_Report_1 <- annotatePeakInBatch(DEAnnotation_Report, AnnotationData=annoData)
    DEAnnotation_Report_2 <- addGeneIDs(DEAnnotation_Report_1, annoOrg)

    return(DEAnnotation_Report_2)

  })

  output$dbareport_Annotation_Table <- renderDataTable({

    DEAnnotation_Report_2 <- DE_Annotation()
    DEAnnotation_Report_3 <- as.data.frame(DEAnnotation_Report_2@elementMetadata)
    DEAnnotation_Report_3[c(-1,-4,-5)]

  })

  #################
  #
  # Peak Annotation Table
  #
  #####################################
  output$ChIPpeakAnno_PeakAnnotation_Table <- renderPrint({
    overlaps.anno<- ChIPpeakAnno_PeakAnnotation()
    head(overlaps.anno)
  })

  output$ChIPpeakAnno_PeakAnnotation_Table_1 <- renderPrint({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]
    peakSummary <- peakData()

    Sample_1.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[1]],AnnotationData=annoData)
    head(Sample_1.anno)[,1:5]
  })

  output$ChIPpeakAnno_PeakAnnotation_Table_2 <- renderPrint({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]
    peakSummary <- peakData()

    Sample_2.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[2]],AnnotationData=annoData)
    head(Sample_2.anno)[,1:5]
  })

  output$ChIPpeakAnno_PeakAnnotation_Table_3 <- renderPrint({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]
    peakSummary <- peakData()

    Sample_3.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[3]],AnnotationData=annoData)
    head(Sample_3.anno)[,1:5]
  })

  output$ChIPpeakAnno_PeakAnnotation_Table_4 <- renderPrint({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]
    peakSummary <- peakData()

    Sample_4.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[4]],AnnotationData=annoData)
    head(Sample_4.anno)[,1:5]
  })

  output$ChIPpeakAnno_PeakAnnotation_Table_5 <- renderPrint({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]
    peakSummary <- peakData()

    Sample_5.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[5]],AnnotationData=annoData)
    head(Sample_5.anno)[,1:5]
  })

  ##################



  # Plot Section ------------------------------------------------------------
  #
  # Volcano Plot
  #
  #####################################

  volcano_plot_MatrixCreate_Plot <- reactive({

    dbaObject <- dbaAnalysis()
    tamoxifen.DB <-  dbareport()
    contrastSelect <- contrastSelect()

    volcanoPlotResult <- DE_VolcanoPlotMatrix(dbaObject,contrast=contrastSelect,
                                              method=input$plot_method,th=input$fdr,bUsePval=FALSE,
                                              fold=input$fold_change,facname="",bLabels=FALSE,
                                              maxLabels=50,bSignificant=TRUE, bFlip=FALSE)
    return(volcanoPlotResult)

  })

  output$dba_Volcano <- renderPlot({

    volcanoPlotDataList <- volcano_plot_MatrixCreate_Plot()
    volcano_plot_matrix <- volcanoPlotDataList$volcano_plot_matrix
    plotTitle <- volcanoPlotDataList$plotTitle
    xLabel <- volcanoPlotDataList$xLabel
    yLabel <- volcanoPlotDataList$yLabel

    ggplot(volcano_plot_matrix,aes(x=log2FC,y=log10FDR)) +
      geom_point(aes(col=Legends),size=3,alpha=0.8) +
      scale_color_manual(values=c("#FBAD01","#36BDE9","#FF0000")) +
      labs(title=plotTitle,x=xLabel,y=yLabel)
  })

  #################
  #
  # Bar Plot
  #
  #####################################
  output$dba_BarPlot <- renderPlot({

    dbaObject <- dbaAnalysis()
    contrastSelect <- contrastSelect()

    DE_BarPlot <- DE_BarPlot(dbaObject, contrastSelect, input)
    DE_BarPlot

  })

  #################
  #
  # PCA Plot
  #
  #####################################
  # PCA plots
  output$dba_PCA <- renderPlot({

    dbaObject <- dbaAnalysis()
    dba.plotPCA(dbaObject,label=DBA_CONDITION)

  })

  #################
  #
  # Heatmap Plot
  #
  #####################################

  output$dba_heatmap <- renderPlot({

    dbaObject <- dbaAnalysis()
    contrastSelect <- contrastSelect()
    hmap <- colorRampPalette(c("blue", "black", "red"))(n = 13)

    readscores <- dba.plotHeatmap(dbaObject, contrast=contrastSelect, correlations=FALSE,
                                  ColAttributes=NULL,
                                  scale="row", colScheme = hmap)
  })

  #################
  #
  # Venn Plot
  #
  #####################################

  # Draw Venn Plot
  output$ChIPpeakAnno_Venn <- renderPlot({

    input$vennPlot
    if (input$vennPlot == 0)return()
    isolate({

      peakSummary <- peakData()
      vennPlotOption_UserInput <- overlapSelection_Return()

      withProgress(message = 'Making plot', value = 0.3,{

        # Judge whether user choose right number
        if(is.null(vennPlotOption_UserInput[["userInput"]])){
          return()
        }else{
          overlaps <- generateOverlapObject()
          makeVennDiagram(overlaps, NameOfPeaks = input$overlapSampleSelection,
                          fill=vennPlotOption_Color(),
                          col=vennPlotOption_Color(),
                          cat.cex=1.2
          )
        }

      })
    })
  })

  ##################
  #
  # Distribution Plot:according to genome
  #
  #####################################
  ChIPpeakAnno_DuplicateDistribution_Plot <- reactive({

    input$duplicateDistributionPlot

    if (input$duplicateDistributionPlot == 0)return()

    isolate({

      peakSummary <<- peakData()

      withProgress(message = 'Making plot', value = 0.3,{

        # Judge whether user choose right number
        peaks <- GRangesList(peakSummary$ChIPpeakAnno_Replicate[input$overlapSampleSelection])
        GED_plot <- genomicElementDistribution(peaks,
                                               TxDb = annoTxDb,
                                               promoterRegion=c(upstream=2000, downstream=500),
                                               geneDownstream=c(upstream=0, downstream=2000),
                                               plot=TRUE)

        return(GED_plot)

      })
    })
  })

  output$ChIPpeakAnno_DuplicateDistribution <- renderPlot({

    ChIPpeakAnno_DuplicateDistribution_Plot()

  })

  #################
  #
  # Peak distribution over different genomic features.
  #
  #####################################
  ChIPpeakAnno_OverlapDistribution_Plot <- reactive({
    input$overlapDistributionPlot
    if (input$overlapDistributionPlot == 0)return()
    isolate({
      overlaps <- generateOverlapObject()
      overlaps_plot <- overlaps$peaklist[[input$overlapSelection]]

      #' check the genomic element distribution for the overlaps
      #' the genomic element distribution will indicates the
      #' the best methods for annotation.
      #' The percentages in the legend show the percentage of peaks in
      #' each category.

      withProgress(message = 'Making plot', value = 0.3,{
        genomicElementDistribution(overlaps_plot,
                                   TxDb = annoTxDb,
                                   promoterRegion=c(upstream=2000, downstream=500),
                                   geneDownstream=c(upstream=0, downstream=2000),
                                   promoterLevel=list(
                                     # from 5' -> 3', fixed precedence 3' -> 5'
                                     breaks = c(-2000, -1000, -500, 0, 500),
                                     labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                "upstream <500b", "TSS - 500b"),
                                     colors = c("#FFE5CC", "#FFCA99",
                                                "#FFAD65", "#FF8E32")))
      })
    })
  })

  output$ChIPpeakAnno_OverlapDistribution <- renderPlot({
    ChIPpeakAnno_OverlapDistribution_Plot()
  })

  #################
  #
  # Pie Plot
  #
  #####################################

  # Overlap Pie plot
  output$ChIPpeakAnno_PeakAnnotation_Pie <- renderPlot({
    overlaps.anno<- ChIPpeakAnno_PeakAnnotation()
    pie1(table(overlaps.anno$insideFeature))
  })

  # Individual Pie plot
  output$ChIPpeakAnno_PeakAnnotation_Pie_1 <- renderPlot({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]
    peakSummary <- peakData()

    Sample_1.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[1]],AnnotationData=annoData)
    pie1(table(Sample_1.anno$insideFeature))
  })

  output$ChIPpeakAnno_PeakAnnotation_Pie_2 <- renderPlot({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]
    peakSummary <- peakData()

    Sample_2.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[2]],AnnotationData=annoData)
    pie1(table(Sample_2.anno$insideFeature))
  })

  output$ChIPpeakAnno_PeakAnnotation_Pie_3 <- renderPlot({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]
    peakSummary <- peakData()

    Sample_3.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[3]],AnnotationData=annoData)
    pie1(table(Sample_3.anno$insideFeature))
  })

  output$ChIPpeakAnno_PeakAnnotation_Pie_4 <- renderPlot({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]
    peakSummary <- peakData()

    Sample_4.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[4]],AnnotationData=annoData)
    pie1(table(Sample_4.anno$insideFeature))
  })

  output$ChIPpeakAnno_PeakAnnotation_Pie_5 <- renderPlot({
    overlapSelection_UserInput <- input$overlapSampleSelection
    overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]
    peakSummary <- peakData()

    Sample_5.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[5]],AnnotationData=annoData)
    pie1(table(Sample_5.anno$insideFeature))
  })

  #################
  #
  # KEGG Pathway Enrichment Plot
  # Subset is the data frame produced after filter, which can fast the react when
  # users change the FDR
  #
  #####################################

  # DE KEGG plot
  output$DE_KEGGplot <- renderPlot({

    KEGG_DE_Pathway <- DE_Annotation_KEGG()
    KEGG_DE_Title <- paste("KEGG pathway enrichment")

    KEGG_DE_Pathway_Subset <- dplyr::filter(KEGG_DE_Pathway, pvalue <= input$fdr)

    KEGG_DE <- KEGGplot(KEGG_DE_Pathway_Subset, KEGG_DE_Title)
    KEGG_DE

  })

  # Overlap KEGG plot
  output$ChIPpeakAnno_PeakAnnotation_KEGGPathwayPlot <- renderPlot({

    KEGG_Overlap_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table()
    KEGG_Overlap_Title <- paste("Overlap KEGG pathway enrichment")

    KEGG_Overlap_Pathway_Subset <- dplyr::filter(KEGG_Overlap_Pathway, pvalue <= input$fdr)

    # Draw the plot
    KEGG_Overlap <- KEGGplot(KEGG_Overlap_Pathway_Subset, KEGG_Overlap_Title)
    KEGG_Overlap

  })

  # Individual KEGG plot
  output$ChIPpeakAnno_PeakAnnotation_PathwayPlot_1 <- renderPlot({

    withProgress(message = 'Making KEGG plot 1', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(1)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_1,"KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)
      KEGG_Individual

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_PathwayPlot_2 <- renderPlot({

    withProgress(message = 'Making KEGG plot 2', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(2)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_2,"KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)
      KEGG_Individual

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_PathwayPlot_3 <- renderPlot({

    withProgress(message = 'Making KEGG plot 3', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(3)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_3,"KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)
      KEGG_Individual

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_PathwayPlot_4 <- renderPlot({

    withProgress(message = 'Making KEGG plot 4', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(4)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_4, "KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)
      KEGG_Individual

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_PathwayPlot_5 <- renderPlot({

    withProgress(message = 'Making KEGG plot 5', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(5)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_5, "KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)
      KEGG_Individual

    })
  })

  #################
  #
  # GO Biological Process Enrichment Plot
  #
  #####################################

  # DE GO Biological Process Plot
  output$DE_GOBiologicalProcessPathwayPlot <- renderPlot({

    GO_DE <- DE_Annotation_GO()
    GObp_DE_Pathway <- GO_DE$bp
    GObp_DE_Title <- paste("GO biological process enrichment")

    GObp_DE_Pathway_Subset <- dplyr::filter(GObp_DE_Pathway, pvalue <= input$fdr)

    GObp_DE_plot <- GObpPlot(GObp_DE_Pathway_Subset, GObp_DE_Title)
    GObp_DE_plot

  })

  # Overlap GO Biological Process Plot
  output$ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot <- renderPlot({

    GO_Overlap <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table()
    GObp_Overlap_Pathway <- GO_Overlap$bp
    GObp_Overlap_Title <- paste("Overlap GO biological process enrichment")

    GObp_Overlap_Pathway_Subset <- dplyr::filter(GObp_Overlap_Pathway, pvalue <= input$fdr)

    # Draw the plot
    GObp_Overlap_plot <- GObpPlot(GObp_Overlap_Pathway_Subset, GObp_Overlap_Title)
    GObp_Overlap_plot

  })

  # Individual GO Biological Process Plot
  output$ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_1 <- renderPlot({

    withProgress(message = 'Making GO plot 1', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(1)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_1,"GO biological process enrichment",sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      GObp_Individual_plot

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_2 <- renderPlot({

    withProgress(message = 'Making GO plot 2', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(2)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_2,"GO biological process enrichment",sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      GObp_Individual_plot

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_3 <- renderPlot({

    withProgress(message = 'Making GO plot 3', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(3)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_3,"GO biological process enrichment",sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      GObp_Individual_plot

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_4 <- renderPlot({

    withProgress(message = 'Making GO plot 4', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(4)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_4,"GO biological process enrichment",sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      GObp_Individual_plot

    })
  })

  output$ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_5 <- renderPlot({

    withProgress(message = 'Making GO plot 5', value = 0.3,{

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(5)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_5,"GO biological process enrichment",sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      GObp_Individual_plot

    })
  })

  #################



  # Download Section --------------------------------------------------------
  #
  # Download DE Report Table
  #
  #####################################
  output$downloadReportData <- downloadHandler(

    filename = function() {

      "DE_Report.csv"

    },
    content = function(file) {

      DE_report <- dbareport()
      write.csv(as.data.frame(DE_report), file, row.names = FALSE )

    }
  )

  #################
  #
  # Download DE Annotation Table
  #
  #####################################
  output$downloadDE_AnnotationReport <- downloadHandler(

    filename = function() {

      "DE_AnnotationReport.csv"

    },
    content = function(file) {

      DEAnnotation_Report_2 <- DE_Annotation()
      DEAnnotation_Report_3 <- as.data.frame(DEAnnotation_Report_2@elementMetadata)

      write.csv(DEAnnotation_Report_3, file, row.names = FALSE )

    }
  )
  #################
  #
  # Download Volcano Plot
  #
  #####################################
  output$downloadVolcanoPlot_DE <- downloadHandler(

    filename = function() {
      paste("DE_Volcanoplot.jpeg")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      volcanoPlotDataList <- volcano_plot_MatrixCreate_Plot()
      volcano_plot_matrix <- volcanoPlotDataList$volcano_plot_matrix
      plotTitle <- volcanoPlotDataList$plotTitle
      xLabel <- volcanoPlotDataList$xLabel
      yLabel <- volcanoPlotDataList$yLabel

      p <- ggplot(volcano_plot_matrix,aes(x=log2FC,y=log10FDR)) +
        geom_point(aes(col=Legends),size=3,alpha=0.8) +
        scale_color_manual(values=c("#FBAD01","#36BDE9","#FF0000")) +
        labs(title=plotTitle,x=xLabel,y=yLabel)

      print(p)

      dev.off()

    }

  )
  #################
  #
  # Download PCA Plot
  #
  #####################################
  output$downloadPCAPlot_DE <- downloadHandler(

    filename = function() {
      paste("DE_PCAplot.jpeg")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      dbaObject <- dbaAnalysis()
      dba.plotPCA(dbaObject,label=DBA_CONDITION)

      dev.off()

    }

  )
  #################
  #
  # Download Heatmap Plot
  #
  #####################################
  output$downloadHeatmapPlot_DE <- downloadHandler(
    filename = function() {

      paste("DE_HeatmapPlot.jpeg")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      dbaObject <- dbaAnalysis()
      contrastSelect <- contrastSelect()
      hmap <- colorRampPalette(c("blue", "black", "red"))(n = 13)

      readscores <- dba.plotHeatmap(dbaObject, contrast=contrastSelect, correlations=FALSE,
                                    ColAttributes=NULL,
                                    scale="row", colScheme = hmap)

      dev.off()

    }
  )

  #################
  #
  # Download Venn Plot
  #
  #####################################
  output$downloadVennPlot_AN <- downloadHandler(

    filename = function() {

      paste("AN_VennPlot.jpeg")

    },
    content = function(file) {

      jpeg(file,width=2500,height=2100,res=300)

      input$vennPlot
      peakSummary <- peakData()
      vennPlotOption_UserInput <- overlapSelection_Return()

      if(is.null(vennPlotOption_UserInput[["userInput"]])){

        return()

      }else{

        overlaps <- generateOverlapObject()
        makeVennDiagram(overlaps, NameOfPeaks = input$overlapSampleSelection,
                        fill=vennPlotOption_Color(),
                        col=vennPlotOption_Color(),
                        cat.cex=1.2
        )}

      dev.off()
    }
  )
  #################
  #
  # Download Genomic Element Distribution Plot
  #
  #####################################

  # Download GED_D plot
  output$downloadDuplicateDistribution_AN <- downloadHandler(
    filename = function() {
      paste("AN_GEDD_Plot.jpeg")
    },
    content = function(file) {
      jpeg(file,width=3200,height=2100,res=300)
      p <- ChIPpeakAnno_DuplicateDistribution_Plot()
      print(p)
      dev.off()
    }
  )

  # Download GED_O plot
  output$downloadOverlapDistributionPlot_AN <- downloadHandler(
    filename = function() {
      paste("AN_GEDO_Plot.jpeg")
    },
    content = function(file) {
      jpeg(file,width=3200,height=2100,res=300)
      p <- ChIPpeakAnno_OverlapDistribution_Plot()
      print(p)
      dev.off()
    }
  )

  #################
  #
  # Download Annotation Table
  #
  #####################################
  output$downloadAnnotatePeaksData <- downloadHandler(
    filename = function() {
      "Overlap_anno.csv"
    },
    content = function(file) {
      overlaps.anno<- ChIPpeakAnno_PeakAnnotation()
      write.csv(as.data.frame(unname(overlaps.anno)), file, row.names = FALSE)
    }
  )

  output$downloadAnnotatePeaksData_1 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]
      paste(overlapSelection_UserInput_1, "anno.csv", sep="_")
    },
    content = function(file) {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]
      peakSummary <- peakData()

      Sample_1.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[1]],AnnotationData=annoData)
      write.csv(as.data.frame(unname(Sample_1.anno)), file, row.names = FALSE)
    }
  )

  output$downloadAnnotatePeaksData_2 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]
      paste(overlapSelection_UserInput_2, "anno.csv", sep="_")
    },
    content = function(file) {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]
      peakSummary <- peakData()

      Sample_2.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[2]],AnnotationData=annoData)
      write.csv(as.data.frame(unname(Sample_2.anno)), file, row.names = FALSE)
    }
  )

  output$downloadAnnotatePeaksData_3 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]
      paste(overlapSelection_UserInput_3, "anno.csv", sep="_")
    },
    content = function(file) {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]
      peakSummary <- peakData()

      Sample_3.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[3]],AnnotationData=annoData)
      write.csv(as.data.frame(unname(Sample_3.anno)), file, row.names = FALSE)
    }
  )

  output$downloadAnnotatePeaksData_4 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]
      paste(overlapSelection_UserInput_4, "anno.csv", sep="_")
    },
    content = function(file) {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]
      peakSummary <- peakData()

      Sample_4.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[4]],AnnotationData=annoData)
      write.csv(as.data.frame(unname(Sample_4.anno)), file, row.names = FALSE)
    }
  )

  output$downloadAnnotatePeaksData_5 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]
      paste(overlapSelection_UserInput_5, "anno.csv", sep="_")
    },
    content = function(file) {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]
      peakSummary <- peakData()

      Sample_5.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[5]],AnnotationData=annoData)
      write.csv(as.data.frame(unname(Sample_5.anno)), file, row.names = FALSE)
    }
  )

  #################
  #
  # Download Pie plod
  #
  #####################################
  output$downloadAnnotatePiePlot <- downloadHandler(
    filename = function() {
      "OverlapPieplot.jpeg"
    },
    content = function(file) {
      jpeg(file,width=2200,height=2000,res=300)
      overlaps.anno<- ChIPpeakAnno_PeakAnnotation()
      pie1(table(overlaps.anno$insideFeature))
      dev.off()
    }
  )

  output$downloadAnnotatePiePlot_1 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]
      paste(overlapSelection_UserInput_1, "Pieplot.jpeg", sep="_")
    },
    content = function(file) {
      jpeg(file,width=2200,height=2000,res=300)
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]
      peakSummary <- peakData()

      Sample_1.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[1]],AnnotationData=annoData)
      pie1(table(Sample_1.anno$insideFeature))
      dev.off()
    }
  )

  output$downloadAnnotatePiePlot_2 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]
      paste(overlapSelection_UserInput_2, "Pieplot.jpeg", sep="_")
    },
    content = function(file) {
      jpeg(file,width=2200,height=2000,res=300)
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]
      peakSummary <- peakData()

      Sample_2.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[2]],AnnotationData=annoData)
      pie1(table(Sample_2.anno$insideFeature))
      dev.off()
    }
  )

  output$downloadAnnotatePiePlot_3 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]
      paste(overlapSelection_UserInput_3, "Pieplot.jpeg", sep="_")
    },
    content = function(file) {
      jpeg(file,width=2200,height=2000,res=300)
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]
      peakSummary <- peakData()

      Sample_3.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[3]],AnnotationData=annoData)
      pie1(table(Sample_3.anno$insideFeature))
      dev.off()
    }
  )

  output$downloadAnnotatePiePlot_4 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]
      paste(overlapSelection_UserInput_4, "Pieplot.jpeg", sep="_")
    },
    content = function(file) {
      jpeg(file,width=2200,height=2000,res=300)
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]
      peakSummary <- peakData()

      Sample_4.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[4]],AnnotationData=annoData)
      pie1(table(Sample_4.anno$insideFeature))
      dev.off()
    }
  )

  output$downloadAnnotatePiePlot_5 <- downloadHandler(
    filename = function() {
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]
      paste(overlapSelection_UserInput_5, "Pieplot.jpeg", sep="_")
    },
    content = function(file) {
      jpeg(file,width=2200,height=2000,res=300)
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]
      peakSummary <- peakData()

      Sample_5.anno<- annotatePeakInBatch(peakSummary$ChIPpeakAnno_Grange[[5]],AnnotationData=annoData)
      pie1(table(Sample_5.anno$insideFeature))
      dev.off()
    }
  )

  #################
  #
  # Download KEGG plot
  #
  #####################################

  # DE KEGG plot
  output$downloadKEGGplot_DE <- downloadHandler(
    filename = function(){

      paste("DE KEGG Enrichment.jpeg")

    },
    content = function(file){

      jpeg(file,width=3200,height=1000,res=300)

      KEGG_DE_Pathway <- DE_Annotation_KEGG()
      KEGG_DE_Title <- paste("KEGG pathway enrichment")

      KEGG_DE_Pathway_Subset <- dplyr::filter(KEGG_DE_Pathway, pvalue <= input$fdr)

      KEGG_DE <- KEGGplot(KEGG_DE_Pathway_Subset, KEGG_DE_Title)
      print(KEGG_DE)

      dev.off()

    }
  )

  # Overlap KEGG plot
  output$downloadKEGGplot <- downloadHandler(
    filename = function() {

      paste("AN_KEGGPlot.jpeg")

    },
    content = function(file) {

      jpeg(file,width=3200,height=1000,res=300)

      KEGG_Overlap_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table()
      KEGG_Overlap_Title <- paste("Overlap KEGG pathway enrichment")

      KEGG_Overlap_Pathway_Subset <- dplyr::filter(KEGG_Overlap_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Overlap <- KEGGplot(KEGG_Overlap_Pathway_Subset, KEGG_Overlap_Title)
      print(KEGG_Overlap)

      dev.off()

    }
  )

  # Individual KEGG plot
  output$downloadKEGGplot_1 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]

      paste(overlapSelection_UserInput_1,"AN_GOPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(1)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_1, "KEGG pathway enrichment", sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)

      print(KEGG_Individual)
      dev.off()

    }
  )

  output$downloadKEGGplot_2 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]

      paste(overlapSelection_UserInput_2,"AN_GOPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(2)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_2, "KEGG pathway enrichment", sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)

      print(KEGG_Individual)

      dev.off()

    }
  )

  output$downloadKEGGplot_3 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]

      paste(overlapSelection_UserInput_3,"AN_GOPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)
      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(3)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_3,"KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)

      print(KEGG_Individual)

      dev.off()

      print(p)
      dev.off()
    }
  )

  output$downloadKEGGplot_4 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]

      paste(overlapSelection_UserInput_4,"AN_GOPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(4)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_4,"KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)

      print(KEGG_Individual)

      dev.off()

    }
  )

  output$downloadKEGGplot_5 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]

      paste(overlapSelection_UserInput_5,"AN_GOPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]

      KEGG_Individual_Pathway <- ChIPpeakAnno_PeakAnnotation_PathwayPlot_Table_i(5)
      KEGG_Individual_Title <- paste(overlapSelection_UserInput_5,"KEGG pathway enrichment",sep=" ")

      KEGG_Individual_Pathway_Subset <- dplyr::filter(KEGG_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      KEGG_Individual <- KEGGplot(KEGG_Individual_Pathway_Subset, KEGG_Individual_Title)

      print(KEGG_Individual)

      dev.off()

    }
  )

  #################
  #
  # GO Biological Process Enrichment Plot
  #
  #####################################

  # DE KEGG plot
  output$downloadGOplot_DE <- downloadHandler(
    filename = function(){

      paste("DE GO biological process Enrichment.jpeg")

    },
    content = function(file){

      jpeg(file,width=3200,height=1000,res=300)

      GO_DE <- DE_Annotation_GO()
      GObp_DE_Pathway <- GO_DE$bp
      GObp_DE_Title <- paste("GO biological process enrichment")

      GObp_DE_Pathway_Subset <- dplyr::filter(GObp_DE_Pathway, pvalue <= input$fdr)

      GObp_DE_plot <- GObpPlot(GObp_DE_Pathway_Subset, GObp_DE_Title)
      print(GObp_DE_plot)

      dev.off()

    }
  )

  # Overlap KEGG plot
  output$downloadGOBiologicalProcessplot <- downloadHandler(
    filename = function() {

      paste("AN_GObpPlot.jpeg")

    },
    content = function(file) {

      jpeg(file,width=3200,height=1200,res=300)

      GO_Overlap <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table()
      GObp_Overlap_Pathway <- GO_Overlap$bp
      GObp_Overlap_Title <- paste("Overlap GO biological process enrichment")

      GObp_Overlap_Pathway_Subset <- dplyr::filter(GObp_Overlap_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Overlap_plot <- GObpPlot(GObp_Overlap_Pathway_Subset, GObp_Overlap_Title)
      print(GObp_Overlap_plot)

      dev.off()

    }
  )

  # Individual GO plot
  output$downloadGOBiologicalProcessplot_1 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]

      paste(overlapSelection_UserInput_1,"AN_GObpPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_1 <- overlapSelection_UserInput[1]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(1)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_1,"GO biological process enrichment",sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      print(GObp_Individual_plot)

      dev.off()

    }
  )

  output$downloadGOBiologicalProcessplot_2 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]

      paste(overlapSelection_UserInput_2,"AN_GObpPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_2 <- overlapSelection_UserInput[2]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(2)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_2,"GO biological process enrichment",sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      print(GObp_Individual_plot)

      dev.off()

    }
  )

  output$downloadGOBiologicalProcessplot_3 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]

      paste(overlapSelection_UserInput_3,"AN_GObpPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_3 <- overlapSelection_UserInput[3]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(3)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_3, "GO biological process enrichment", sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      print(GObp_Individual_plot)

      dev.off()

    }
  )

  output$downloadGOBiologicalProcessplot_4 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]

      paste(overlapSelection_UserInput_4, "AN_GObpPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_4 <- overlapSelection_UserInput[4]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(4)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_4, "GO biological process enrichment", sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      print(GObp_Individual_plot)

      dev.off()

    }
  )

  output$downloadGOBiologicalProcessplott_5 <- downloadHandler(
    filename = function() {

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]

      paste(overlapSelection_UserInput_5, "AN_GObpPlot.jpeg", sep="_")

    },
    content = function(file) {

      jpeg(file,width=3200,height=2100,res=300)

      overlapSelection_UserInput <- input$overlapSampleSelection
      overlapSelection_UserInput_5 <- overlapSelection_UserInput[5]

      GO_Individual <- ChIPpeakAnno_PeakAnnotation_GOBiologicalProcessPathwayPlot_Table_i(5)
      GObp_Individual_Pathway <- GO_Individual$bp
      GObp_Individual_Title <- paste(overlapSelection_UserInput_5, "GO biological process enrichment", sep=" ")

      GObp_Individual_Pathway_Subset <- dplyr::filter(GObp_Individual_Pathway, pvalue <= input$fdr)

      # Draw the plot
      GObp_Individual_plot <- GObpPlot(GObp_Individual_Pathway_Subset, GObp_Individual_Title)
      print(GObp_Individual_plot)

      dev.off()

    }
  )

  #################

}



