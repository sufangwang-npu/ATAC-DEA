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

  If you want to use this app, please cite as: ATAC-DEA: a web-based ATAC-seq data differential expression and annotation analysis application (in preparation)



## Instruction

### Get a quick start of ATAC-DEA

ATAC-DEA is a web-based platform to help you analyze ATAC-seq data and plot high quality figures. The ATAC-DEA allows users to visualize differential expression(DE) peaks and do annotation. Explore the app's features with the example data set pre-loaded. Upload your genes Expression data first,then submit your data.



### Data Requirements

![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/CountR.jpg)


1. DEList[.csv]:sample sheet contains the main information of each sample
2. mapping results files[.bam]
3. peak calling results files[.bed]



### Example DEList format

![image](https://github.com/sufangwang-npu/ATAC-DEA/blob/main/WWW/DEList.jpg)

**SampleID:** ID of each sample

**Tissue/Factor/Condition/Replicate:** The variables of each sample

**bamReads:** Path of mapping results files[.bam] (==necessary==)

**ControlID:** ID of the control data[.bam] (not necessary)

**bamControl:** Path of the control data[.bam] (not necessary)

**Peaks:** Path of peak calling results files[.bed] (necessary)

**PeakCaller:** Suffix of peak calling results files[.bed] (necessary)



### Data Pretreatment

The three files need to be processed in order to get the right format which can be operated by ATAC-DEA. We wrote a script, called readCount.R, for user to easily do pretreatment. It will take the three needed files and output a peak data collection file which is used as the input of ATAC-DEA.

1. change the directory to your own
2. choose the count option according to your time and computer
3. save the peak collection result file to your working directory



