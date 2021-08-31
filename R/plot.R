library(dplyr)
library(ggplot2)

#' DE_VolcanoPlotMatrix
#'
#'
#' Process the dbaobject data to plot volcano
#'
#' @param pv, DBAobject
#' @param contrast, integer, contrast users selected
#' @param method, 'DESeq2', 'EdgeR', 'Both'
#' @return A list contains volcano plot matrix, plot title, x-label and y-label
#' @export
#'
#' @examples
#'     volcanoPlotResult <- DE_VolcanoPlotMatrix(dbaObject,contrast=contrastSelect,
#'                                               method=input$plot_method,th=input$fdr,bUsePval=FALSE,
#'                                               fold=input$fold_change,facname="",bLabels=FALSE,
#'                                               maxLabels=50,bSignificant=TRUE, bFlip=FALSE)

DE_VolcanoPlotMatrix <- function(pv,contrast,method='DESeq2', th=0.05,
                                 bUsePval=FALSE,fold=0,facname="",
                                 bLabels=FALSE,maxLabels=50,
                                 bSignificant=TRUE, bFlip=FALSE){
  Legend <- Fold <- NULL

  if(missing(contrast)){
    contrast <- 1:length(pv$contrasts)
  } else {
    if(contrast > length(pv$contrasts)) {
      stop('Specified contrast number is greater than number of contrasts',call.=FALSE)
      return(NULL)
    }
  }

  for(con in 1:length(contrast)) {
    conrec <- pv$contrasts[[contrast[con]]]
    name1 <- conrec$name1
    name2 <- conrec$name2
    if(bFlip) {
      name1 <- conrec$name2
      name2 <- conrec$name1
    }
    for(meth in method) {
      res <- dba.report(pv,contrast,bUsePval=TRUE,th=100,bNormalized=TRUE,precision=0)
      res <- DE_VolcanoPlotMatrix_ExcludeNonAnnoPeaks(res)

      if(!is.null(res)) {
        if(bUsePval) {
          vals <- res$"p-value"
          idx  <- vals <= th
          tstr <- "p"
          res<- as.data.frame(res)
          res = mutate(res,
                       Legend=ifelse(res$"p-value"<=th,
                                     sprintf(" p-val<=%1.4f",th),
                                     sprintf(" p-val >%1.4f",th)))
        } else {
          vals <- res$FDR
          idx  <- vals <= th
          tstr <- "FDR"
          res<- as.data.frame(res)
          res = mutate(res,
                       Legend=ifelse(res$FDR<th,
                                     sprintf(" FDR<=%1.4f",th),
                                     sprintf(" FDR >%1.4f",th)))
        }

        res$Legend[idx & abs(res$Fold) < fold] <-
          sprintf("abs(log2FC)<%1.2f",2^fold)
        idx <- idx & abs(res$Fold) >= fold

        sigSites <- res[idx,]
        rownames(sigSites) <- 1:sum(idx)

        res <- cbind(0,res)
        colnames(res)[1] <- "SiteNum"
        res[idx,1] <- 1:sum(idx)

        constr <- getContrastString(conrec)
        plotTitle <- sprintf('%s Contrast: %s  [%s %s<=%1.3f',
                             facname, constr,sum(idx),tstr,th)
        if(fold>0) {
          plotTitle <- sprintf("%s & abs(FoldChange)>=%1.2f]",
                               plotTitle, 2^fold)
        } else {
          plotTitle <- sprintf("%s]",plotTitle)
        }
        if(is.null(conrec$name2)) {
          xLabel <- "log Fold Change"
        } else {
          xLabel <- sprintf('log Fold Change [log2(%s) - log2(%s)]',name1,name2)
        }

        yLabel <- sprintf("-log10(%s)",tstr)

        volcano_plot_matrix <- data.frame(log2FC=res$Fold, log10FDR=-log10(vals),Legend=res$Legend)

        Legends <- c()
        for(i in 1:length(volcano_plot_matrix[[1]])){
          if(volcano_plot_matrix[i,1]>=0 & identical(volcano_plot_matrix[i,3], sprintf(" FDR<=%1.4f",th))){
            Legends <- append(Legends, "Up-regulate")
          }
          else if(volcano_plot_matrix[i,1]<=0 & identical(volcano_plot_matrix[i,3], sprintf(" FDR<=%1.4f",th))){
            Legends <- append(Legends, "Down-regulate")
          }
          else{
            Legends <- append(Legends, sprintf(" FDR>=%1.4f",th))
          }
        }


        volcano_plot_matrix <- cbind(volcano_plot_matrix, Legends)


      }
    }
  }

  volcanoPlotDataList <- list(volcano_plot_matrix=volcano_plot_matrix,plotTitle=plotTitle,xLabel=xLabel, yLabel=yLabel)

  return(volcanoPlotDataList)
}


#' getContrastString
#'
#'
#' Get the Name of contrast in the plot
#'
#' @param conrec
#' @return Name of contrast
#' @export
#'
#' @examples
#'     constr <- getContrastString(conrec)

getContrastString <- function(conrec) {
  if(!is.null(conrec$contrastType)) {
    constr <- conrec$contrast
    if(conrec$contrastType=='bycolumn') {
      constr <- paste(constr,collapse=":")
    } else if(conrec$contrastType=='byresults1') {
      constr <- constr[[1]]
    } else if(conrec$contrastType=='byresults2') {
      constr <- paste(constr[[1]],"vs.",constr[[2]])
    } else if(length(constr) == 2) {
      constr <- paste(constr[1],"vs.",constr[2])
    } else if(length(constr) == 3)  {
      constr <- paste(constr[2],"vs.",constr[3])
    }
  } else {
    constr <- paste(conrec$name1, "vs.", conrec$name2)
  }
  return(constr)
}


#' DE_VolcanoPlotMatrix_ExcludeNonAnnoPeaks
#'
#'
#' Remove peaks that cannot be annotated
#'
#' @param volcanoPlotReport, GRanges data contains the result from dba.report
#' @return DE peaks after exclusion
#' @export
#'
#' @examples
#'     res <- DE_VolcanoPlotMatrix_ExcludeNonAnnoPeaks(res)

DE_VolcanoPlotMatrix_ExcludeNonAnnoPeaks <- function(volcanoPlotReport){

  DEAnnotation_Report <- as.data.frame(volcanoPlotReport)

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

    DE_ExcludeNonAnnoPeaks_Report <- DEAnnotation_Report

  }else{

    DE_ExcludeNonAnnoPeaks_Report <- data.frame(Chr=firstLine,DEAnnotation_Report[-exclude,][-1])

  }


  return(DE_ExcludeNonAnnoPeaks_Report)

}


#' DE_BarPlot
#'
#'
#' plot bar plot of up/down-regulated peaks
#'
#' @param dbaObject, DBAobject
#' @param contrastSelect, integer, contrast users selected
#' @param input, input of shiny
#' @return bar plot
#' @export
#'
#' @examples
#'     DE_BarPlot <- DE_BarPlot(dbaObject, contrastSelect, input)

DE_BarPlot <- function(dbaObject, contrastSelect, input){

  res <- dba.report(dbaObject,contrast=contrastSelect,bUsePval=TRUE,th=100,bNormalized=TRUE,precision=0)
  vals_Fold <- res$Fold
  vals_FDR <- res$FDR
  th_fold_change <- input$fold_change
  th_FDR <- input$fdr
  idx_Down  <- vals_Fold <= -th_fold_change & vals_FDR <= th_FDR
  idx_Up <- vals_Fold >= th_fold_change & vals_FDR <= th_FDR

  upR <- length(which(idx_Up==1))
  downR <- length(which(idx_Down==1))

  bindingSite <- data.frame(Peak=c("Up-regulation", "down-regulation"),
                            Num=c(upR, downR))
  order <- sort(bindingSite$Peak, index.return=TRUE, decreasing = TRUE)
  bindingSite$Peak <- factor(bindingSite$Peak, levels = bindingSite$Peak[order$ix])

  barplot_color <- c("#FF0000", "#36BDE9")

  p <- ggplot(data=bindingSite, aes(Peak, Num, fill=Peak))+
    geom_bar(stat="identity", width=0.8, colour="black", size=0.25, alpha=1)+
    geom_text(aes(x=Peak,y=Num+(upR+downR)/20,label=Num),color="black", position=position_dodge(width=0.9),
              show.legend = F) +
    theme(panel.grid.major = element_line(linetype = "blank"), legend.position="none")+
    scale_fill_manual(values=barplot_color)

  return(p)

}


#' KEGGplot
#'
#'
#' plot KEGG enrichment pathway plot
#'
#' @param pathway, KEGG info-table from getEnrichedPATH()
#' @param KEGGtitle, plot title
#' @return KEGG plot
#' @export
#'
#' @examples
#'     KEGG_DE <- KEGGplot(KEGG_DE_Pathway, KEGG_DE_Title)

KEGGplot <- function(pathway, KEGGtitle){

  order <-sort(pathway$count.InDataset, index.return=TRUE, decreasing = TRUE)
  pathway$path.term <- factor(pathway$path.term, levels = unique(pathway$path.term[order$ix]))

  p <- ggplot(pathway,aes(count.InDataset,path.term))+
    geom_point()+
    geom_point(aes(size=count.InDataset))+
    geom_point(aes(size=count.InDataset,color=-1*log10(pvalue)))+
    scale_colour_gradient(low="blue",high="red")+
    labs(color = expression(-log[10](pvalue)),size = "Gene number",x = "Rich factor",y = "Pathway name",
         title = KEGGtitle)

  return(p)

}


#' GObpPlot
#'
#'
#' plot GO biological process plot
#'
#' @param pathway, KEGG info-table from getEnrichedPATH()
#' @param GOtitle, plot title
#' @return GO biological process plot
#' @export
#'
#' @examples
#'     GObp_Individual_plot <- GOplot(GObp_Individual_Pathway, GObp_Individual_Title)

GObpPlot <- function(pathway, GOtitle){

  order <- sort(pathway$count.InDataset, index.return=TRUE, decreasing = TRUE)
  pathway$go.term <- factor(pathway$go.term, levels = pathway$go.term[order$ix])

  p <- ggplot(pathway,aes(count.InDataset, go.term))+
    geom_point()+
    geom_point(aes(size=count.InDataset))+
    geom_point(aes(size=count.InDataset,color=-1*log10(pvalue)))+
    scale_colour_gradient(low="blue",high="red")+
    labs(color = expression(-log[10](pvalue)),size = "Gene number",x = "Rich factor",y = "Pathway name",
         title = GOtitle)

  return(p)

}

