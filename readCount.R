# readCount.R
# users can use this file to create peak data collection file

library(DiffBind)

# the path of your DEList and the directory
DEList_path = ".../ATAC-seqData/DEList.csv"
working_directory = ".../ATAC-seqData"

DEList_path = "F:/Experimental Data/ATAC-seq/ATAC-seq Data/DEList.csv"
working_directory = "F:/Experimental Data/ATAC-seq/ATAC-seq Data"

setwd(working_directory)
samples <- read.csv(DEList_path)

DEList <- dba(sampleSheet=samples)

# This option is slower but uses the more standard counting function. If you are
# running this app on supercomputer, we recommend this option
DEList <- dba.count(DEList, bUseSummarizeOverlaps = TRUE)

# This option is faster but uses the less standard counting function. If you are
# running this app on PC, we recommend this option
DEList <- dba.count(DEList, bUseSummarizeOverlaps = FALSE)

# save as peak data collection file
save(DEList,file="peak_data_collection")
