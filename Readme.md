# ATAC-DEA
Differential Expression and Annotation of peak in ATAC-seq analyses application (ATAC-DEA)
Copyright (C) 2021 sufangwang-npu  



## Information for the ATAC-DEA App
ATAC-DEA is a web-based platform to help you analyze ATAC-seq data and plot high quality figures. The ATAC-DEA allows users to visualize differential expression(DE) peaks and do annotation. Explore the app's features with the example data set pre-loaded. Upload your genes Expression data first,then submit your data. To use ATAC-DEA, users could do in two ways: (1) run ATAC-DEA directly through the webpage; (2) download the source code and drive the package in R or RStudio locally. 

ATAC-DEA website: http://www.atac-dea.xyz:3838/ATAC-DEA

ATAC-DEA source code: https://github.com/sufangwang-npu/ATAC-DEA


Operating system(s): Platform independent

Programming language: R(>= 4.0)


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

  If you want to use this app, please cite as: Shilong Zhang, Sufang Wang. ATAC-DEA: a web-based ATAC-seq data differential expression and annotation analysis application (under revision)



## Instruction


### Get a quick start of ATAC-DEA

  ATAC-DEA is a web-based platform to help you analyze ATAC-seq data and plot high quality figures. The ATAC-DEA allows users to visualize differential expression(DE) peaks and do annotation. Explore the app's features with the example data set pre-loaded. Upload your genes Expression data first,then submit your data.  


### Required files in preparation
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


### ATAC-DEA Data Preparation Guide
#### 1.	Organize your data files
To do data preparation, three files are needed. We recommend organizing all of the data files in one directory, and placing mapping results files (generally .bam file) in ‘./reads’; placing peak calling results files (generally .bed file) in ‘./peaks’. We use some example files in ‘ATAC-DEA/preparation/extra’ to hint.

#### 2. Install dependency packages
Before the work, some dependency packages needs to be installed
```
install.packages(c("shiny", "shinyFiles", "stringr"))
                      
if (!require("BiocManager", quietly = TRUE))
                      
install.packages("BiocManager")
                      
BiocManager::install("DiffBind")
```

#### 3.	Run preparation shiny application
We design a local shiny app to help you do preparation. Firstly, you can download the source of ATAC-DEA including Data Preparation Application in Github ATAC-DEA (extra example preparation files: http://59.110.11.223:3838/ATAC-DEA_DataPretreatment.zip). After that, you can run these commands in terminal to start the app:
```
#cd ./ATAC-DEA/DataPreparation
                      
#R -e 'shiny::runApp ()'
```
After the all starts successfully, type the url (http://127.0.0.1:xxxx) in your browser.

Or, you can run this app in Rstudio.
#### 4.	The panel of DataPreparation
DataPreparation helps you build up the sample sheet that DiffBind needs.

1)	Open ‘Create DEList’ tab. If you have completed DEList, you can directly upload it and you can also edit it in this panel.

2)	Set the work space, if the datapath in your DEList is relative path. The final result file ‘peak_data_collection’ file will also be sived in the work space

3)	Input the information of your sample. The DEList in ATAC-DEA/DataPreparation/extra shows what each variable represents.

4)	Select the mapping results file and peak calling results file of each sample. These two buttons only record the relative paths, so don't need to worry about the size of files.

5)	Select ‘Insert’ to add a new line and select ‘remove’ to remove the latest line.

6)	You can check your DEList here

7)	When both DEList and work space are done, you can select ‘complete’ to go to the next step. If one of these two steps has not been done, ‘complete’ button will not response.
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/1.jpg)

8) The page will automatically redirect to readCount tab. You can check your DEList again. If there is no mistake, you can press ‘Do Analysis!’ to start counting the reads. It will take some time, depending on your computer.
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/2.png)

9) When counting is done, a message window will pop out and a file named peak_data_collection will be saved in your work space. Now you can upload this file in ATAC-DEA to explore your data.
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/3.png)

### ATAC-DEA Quick-start Tutorial
This quick-start tutorial will guide you through using the example data to do DE & annotation analysis, and visualize the result.

Open http://www.atac-dea.xyz:3838/ATAC-DEA

1.	In the ‘Upload Your Data’ tab, you can select to use example or your own data

1)	Example data: no files need to be uploaded. ATAC-DEA will use data of breast cancer cell to help you explore it.

2)	Your own data: upload your peak collection result file generated in preparation step and peak calling results files of your data. You can check your data in ‘Check your Sample Sheet’ box.

2.	Select the organism of your data: ATAC-DEA now provides ‘Homo sapiens’ and ‘Mus musculus’ to select. When all above is done, select ‘Go Next’ button
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/Tutorial1.jpg)

3.	A new sub tab ‘Contrast Design’ will be created and you can design your contrast model. You can use default option to explore ATAC-DEA. The detail of other options is listed in the paper. After that, select ‘Do analysis’ to process the analysis and it will take some time.

4.	ATAC-DEA will analysis all possible contrast results in your contrast model and show them at the bottom of the box for you to select. Because the table is interactive, you can directly select one for further analysis.
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/Tutorial2.jpg)

5.	Four new tabs will be shown in sidebar, and you can explore the analysis results by selecting ‘Data Report’, ‘DE Analysis’ and ‘Peak Annotation’
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/Tutorial3.jpg)

6.	Filter is used to set threshold and method (DESeq2 and edgeR are provided) in volcano plot.

7.	Volcano plot is interactive, so you can use brush to select some dots and check their information.
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/Tutorial4.jpg)

8.	Explore the result of DE analysis in ‘Volcano’, ‘PCA’, ‘Heatmap’ and ‘KEGG&GO’ sub tabs
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/Tutorial5.jpg)

9.	In ‘Peak Annotation’ tab, you should select the samples you would like to annotate (more than 2 and no more than 5 otherwise ATAC-DEA cannot get overlaps)

10.	ATAC-DEA will show all overlaps according to the samples you selected. After that select ‘Do Analysis!’ and select the sub tabs to explore the results of annotation.
![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/Tutorial6.jpg)

  

### Data Input
  ATAC-DEA takes two files as input: Users can start analyses with uploading **the peak data collection file (generated from data pretreatment)** and **peak calling results files(.bed)** in “Data Upload” panel of ATAC-DEA. As users finish uploading, DEList will be displayed at the bottom of the page for users to check if the data is completed and correct.   

  

### Data Analysis
  ATAC-DEA uses DiffBind to perform differential expression (DE) analysis. The first step of DE analysis is establishment of contrast model. ATAC-DEA provides two ways: (1) set up a specific contrast model according to the factor selected by users. Single factor and multiple factors analysis are both allowed in ATAC-DEA, and it will establish the contrast model on the grounds of user’s choice; (2) set up all-possible contrasts, implementing by ATAC-DEA automatically. ATAC-DEA will search all possible contrasts, arrange them in an interactive table and users can select one of them to establish the contrast model.
  After the establishment of contrast model, users can press “Do analysis” button and explore the analysis results by selecting the tabs including Data Report, DE Analysis and Peak Annotation in the sidebar of ATAC-DEA. In DE analysis, ATAC-DEA can find loci exhibiting significant differences between different treatment conditions by using DESeq2 or edgeR. 
  ChIPpeakAnno is used to process peak annotation. ATAC-DEA enables the binding sites annotation from human and mouse, which can show the relevant genes of each peak and output figures to illustrate the distribution of genomic elements. For each replicate, the significance of the overlap genes is determined and ATAC-DEA will merge peaks across replicates to obtain the distances between peak location and nearest transcription start site (TSS) and nearest genes, based on their genomic location (i.e., intron, exon, promoter, untranslated regions (UTRs)). Existing annotation packages, such as GenomicFeatures and BSgenome, can be used to annotate the peak set according to the result of merging, and by using ChIPpeakAnno’s getEnrichedGO and getEnrichedPATH functions, analysis for over-represented gene ontology terms and KEGG pathways can be accomplished, respectively.  

  

### Data Download
  All tables (DE report table and Annotation report table) and figures created from ATAC-DEA can be directly downloaded from the application. It supports .csv formats for tables, and .jpeg formats for figures with high qualities.  