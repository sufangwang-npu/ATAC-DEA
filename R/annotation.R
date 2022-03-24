#' addChr2DEReport
#'
#'
#' If the raw data use '1', '2' ... as chromosome, which cannot be recognized by EnsDb,
#' addChr2DEReport can convert them into 'chr1','chr2' ...
#'
#' @param DEAnnotation_Report_noAnno, GRanges data contains the result from dba.report
#' @return Differential expression peaks in GRanges type.
#' @export
#'
#' @examples
#'     DEAnnotation_Report <- addChr2DEReport(DEAnnotation_Report_noAnno)

addChr2DEReport <- function(DEAnnotation_Report_noAnno){

  DEAnnotation_Report <- as.data.frame(DEAnnotation_Report_noAnno)

  chromosome <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X","Y","M")
  firstLine <- c()
  exclude <- c()
  for(i in DEAnnotation_Report[[1]]){
    if(i %in% chromosome){
      j <- paste("chr", i, sep="")
      firstLine <- append(firstLine, j)
    }else if(i == "MT"){
      j <- "chrM"
      firstLine <- append(firstLine, j)
    }else if(i %in% annoData@seqnames@values){
      j <- i
      firstLine <- append(firstLine, j)
    }else{
      j <- i
      exclude_num <- which(DEAnnotation_Report[[1]] == j)
      exclude <- append(exclude, exclude_num)
      exclude <- unique(exclude)
    }

  }

  peakNameLine <-c()
  for(i in 1:length(DEAnnotation_Report[[1]])){
    j <- paste("DE_peak_", ROWNAMES(DEAnnotation_Report)[i], sep="")
    peakNameLine <- append(peakNameLine, j)
  }

  if(is.null(exclude)){

    DE_peak <- data.frame(DEAnnotation_Report[,1:3],Name=peakNameLine,DEAnnotation_Report[,4])

  }else{

    DE_peak <- data.frame(Chr=firstLine, Start=DEAnnotation_Report[[2]][-exclude],
                          End=DEAnnotation_Report[[3]][-exclude],Name=peakNameLine[-exclude],
                          Score=DEAnnotation_Report[[4]][-exclude])

  }


  DE_peak_GRangs <- toGRanges(DE_peak, format="BED", header=FALSE)

  return(DE_peak_GRangs)

}



