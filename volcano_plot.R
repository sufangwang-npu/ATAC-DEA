library(dplyr)
library(ggplot2)

DBAplotVolcano <- function(pv,contrast,method='DESeq2', th=0.05,
                           bUsePval=FALSE,fold=0,facname="",
                           bLabels=FALSE,maxLabels=50,
                           dotSize=3,bSignificant=TRUE, bFlip=FALSE,
                           xrange,yrange) {
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
          plotTitle <- sprintf("%s & abs(log2FC)>=%1.2f]",
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

        volcano_plot_matrix <- data.frame(threshold=res$Fold, log10Value=-log10(vals),Legend=res$Legend)

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

        volcano_plot_matrix_threshold <- cbind(volcano_plot_matrix, Legends)



        volcanoPlot <- DBAplotVolcano_Plot(volcano_plot_matrix_threshold, dotSize, plotTitle, xLabel,
                       yLabel, bLabels, sigSites, idx, maxLabels, Fold, vals, th)



      }
    }
  }
  #return(sigSites[,-10])
  volcanoPlotResult <- list(volcano_plot_matrix = volcano_plot_matrix, volcanoPlot = volcanoPlot)
  return(volcanoPlotResult)
}

# Draw the volcano plot
DBAplotVolcano_Plot <- function(volcano_plot_matrix_threshold, dotSize, plotTitle, xLabel,
                           yLabel, bLabels, sigSites, idx, maxLabels, Fold, vals, th){


  p <- ggplot(volcano_plot_matrix_threshold,aes(x=threshold,y=log10Value)) +
    geom_point(aes(col=Legends),size=dotSize,alpha=0.8) +
    scale_color_manual(values=c("#FBAD01","#36BDE9","#FF0000")) +
    labs(title=plotTitle,x=xLabel,y=yLabel)

  if(bLabels) {
    maxLabels <- min(sum(idx),maxLabels)
    if(maxLabels > 0) {
      xx <-  which(idx)[1:maxLabels]
      p <- p + geom_text_repel(data=sigSites[1:maxLabels,],
                               aes(x=Fold,
                                   y = -log10(vals[xx]),
                                   label=rownames(sigSites)[1:maxLabels]))
    }
  }
  return(p)
}

# Get the Name of contrast in the plot
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
