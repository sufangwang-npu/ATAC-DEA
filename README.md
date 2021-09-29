# ATAC-DEA

Differential Expression and Annotation of peak in ATAC-seq analyses application (ATAC-DEA)
Copyright (C) 2021 sufangwang-npu  



## Information for the ATAC-DEA App

Code canbe found on github:https://github.com/sufangwang-npu/ATAC-DEA.  
To run this app locally on your machine,download R or Rstudio and run the following command once to set up the environment:  
install.packages(c("shiny","shinydashboard","BiocManager","DiffBind","shinyjs","bslib","ChIPpeakAnno","ggplot2","dplyr","DT","reactome.db","RColorBrewer","TxDb.Hsapiens.UCSC.hg19.knownGene","EnsDb.Hsapiens.v75","org.Hs.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene","EnsDb.Mmusculus.v79","org.Mm.eg.db"))  

If you were ready for this packages, You may now run the shiny app with just one command in R:  
library("shiny")  
runApp("ATAC-DEA")  
Or,  
shiny::runGitHub("ATAC-DEA","sufangwang-npu")   
  
  


## Licensing

ATAC-DEA
Differential Expression and Annotation of peak in ATAC-seq analyses application

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  
This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  
You may contact the author of this code, Sufang Wang, at <sufangwang@nwpu.edu.cn>  
  


## Citation

  If you want to use this app, please cite as: Shilong Zhang, Sufang Wang. ATAC-DEA: a web-based ATAC-seq data differential expression and annotation analysis application (in preparation)
  


## Instruction

### Get a quick start of ATAC-DEA

ATAC-DEA is a web-based platform to help you analyze ATAC-seq data and plot high quality figures. The ATAC-DEA allows users to visualize differential expression(DE) peaks and do annotation. Explore the app's features with the example data set pre-loaded. Upload your genes Expression data first,then submit your data.  


### Data Preparation and Pretreatment

#### Data requirements

![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/CountR.jpg)


1. DEList[.csv]:sample sheet contains the main information of each sample
2. mapping results files[.bam]
3. peak calling results files[.bed]



#### Example DEList format

![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/DEList.jpg)

**SampleID:** ID of each sample

**Tissue/Factor/Condition/Replicate:** The variables of each sample

**bamReads:** Path of mapping results files[.bam] (***necessary***)

**ControlID:** ID of the control data[.bam] (not necessary)

**bamControl:** Path of the control data[.bam] (not necessary)

**Peaks:** Path of peak calling results files[.bed] (***necessary***)

**PeakCaller:** Suffix of peak calling results files[.bed] (***necessary***)



#### Process the three files

The three files need to be processed in order to get the right format which can be operated by ATAC-DEA. We wrote a script, called readCount.R, for user to easily do pretreatment. It will take the three needed files and output a peak data collection file which is used as the input of ATAC-DEA.

1. change the directory to your own
2. choose the count option according to your time and computer
3. save the peak collection result file to your working directory  
  
  
  
  
### Data Input
ATAC-DEA takes two files as input: Users can start analyses with uploading **the peak data collection file (generated from data pretreatment)** and **peak calling results files(.bed)** in “Data Upload” panel of ATAC-DEA. As users finish uploading, DEList will be displayed at the bottom of the page for users to check if the data is completed and correct.   

  
  
### Data Analysis
ATAC-DEA uses DiffBind to perform differential expression (DE) analysis. The first step of DE analysis is establishment of contrast model. ATAC-DEA provides two ways: (1) set up a specific contrast model according to the factor selected by users. Single factor and multiple factors analysis are both allowed in ATAC-DEA, and it will establish the contrast model on the grounds of user’s choice; (2) set up all-possible contrasts, implementing by ATAC-DEA automatically. ATAC-DEA will search all possible contrasts, arrange them in an interactive table and users can select one of them to establish the contrast model.
After the establishment of contrast model, users can press “Do analysis” button and explore the analysis results by selecting the tabs including Data Report, DE Analysis and Peak Annotation in the sidebar of ATAC-DEA. In DE analysis, ATAC-DEA can find loci exhibiting significant differences between different treatment conditions by using DESeq2 or edgeR. 
ChIPpeakAnno is used to process peak annotation. ATAC-DEA enables the binding sites annotation from human and mouse, which can show the relevant genes of each peak and output figures to illustrate the distribution of genomic elements. For each replicate, the significance of the overlap genes is determined and ATAC-DEA will merge peaks across replicates to obtain the distances between peak location and nearest transcription start site (TSS) and nearest genes, based on their genomic location (i.e., intron, exon, promoter, untranslated regions (UTRs)). Existing annotation packages, such as GenomicFeatures and BSgenome, can be used to annotate the peak set according to the result of merging, and by using ChIPpeakAnno’s getEnrichedGO and getEnrichedPATH functions, analysis for over-represented gene ontology terms and KEGG pathways can be accomplished, respectively.  

  
  
### Data Download
All tables (DE report table and Annotation report table) and figures created from ATAC-DEA can be directly downloaded from the application. It supports .csv formats for tables, and .jpeg formats for figures with high qualities.  


