library(gplots)
library(RColorBrewer)
library(ggplot2)
library(tibble)
library(scales)

load("~/tyson/VikasDifferentialExpressionControlvs1hrPlasma/Vikas.Rdata")
setwd("~/tyson/VikasDifferentialExpressionControlvs1hrPlasma/")
allDiff <- res.sig.1hr
index=which(rownames(allDiff)=="ENSG00000007952" |
              rownames(allDiff)=="ENSG00000074771" |
              rownames(allDiff)=="ENSG00000086991" |
              rownames(allDiff)=="ENSG00000165168" |
              rownames(allDiff)=="ENSG00000285441" |
              rownames(allDiff)=="ENSG00000143878" |
              rownames(allDiff)=="ENSG00000142168" |
              rownames(allDiff)=="ENSG00000112096" |
              rownames(allDiff)=="ENSG00000109572" |
              rownames(allDiff)=="ENSG00000141837" |
              rownames(allDiff)=="ENSG00000151067" |
              rownames(allDiff)=="ENSG00000285479" |
              rownames(allDiff)=="ENSG00000157388" |
              rownames(allDiff)=="ENSG00000198216" |
              rownames(allDiff)=="ENSG00000006283" |
              rownames(allDiff)=="ENSG00000196557" |
              rownames(allDiff)=="ENSG00000007402" |
              rownames(allDiff)=="ENSG00000165995" |
              rownames(allDiff)=="ENSG00000108878" |
              rownames(allDiff)=="ENSG00000178363" |
              rownames(allDiff)=="ENSG00000178372" |
              rownames(allDiff)=="ENSG00000169885" |
              rownames(allDiff)=="ENSG00000058404" |
              rownames(allDiff)=="ENSG00000240583" |
              rownames(allDiff)=="ENSG00000167580" |
              rownames(allDiff)=="ENSG00000165272" |
              rownames(allDiff)=="ENSG00000171885" |
              rownames(allDiff)=="ENSG00000161798" |
              rownames(allDiff)=="ENSG00000086159" |
              rownames(allDiff)=="ENSG00000165269" |
              rownames(allDiff)=="ENSG00000103375" |
              rownames(allDiff)=="ENSG00000103569" |
              rownames(allDiff)=="ENSG00000143595" |
              rownames(allDiff)=="ENSG00000178301" |
              rownames(allDiff)=="ENSG00000047936" |
              rownames(allDiff)=="ENSG00000136888" |
              rownames(allDiff)=="ENSG00000005381" |
              rownames(allDiff)=="ENSG00000122986" |
              rownames(allDiff)=="ENSG00000128340" |
              rownames(allDiff)=="ENSG00000146648" |
              rownames(allDiff)=="ENSG00000135679" |
              rownames(allDiff)=="ENSG00000198625" |
              rownames(allDiff)=="ENSG00000125995" |
              rownames(allDiff)=="ENSG00000198015" |
              rownames(allDiff)=="ENSG00000114026" |
              rownames(allDiff)=="ENSG00000120129" |
              rownames(allDiff)=="ENSG00000100292" |
              rownames(allDiff)=="ENSG00000111206" )

Genes=allDiff[index, ]
GeneNames=row.names(Genes)
FriendlyName=Genes$symbol
log2Genes=log2.counts[GeneNames, ]
#log2.counts=log2.counts[,1:117]
condition=rep(c("Control"), times=4)
condition=as.character(c(condition, rep(c("Plasma 1 hr"), times=4)))
plotList <- c()
log2.counts=log2.counts[,1:8]
for (x in 1:length(GeneNames)) {
  vals=as.numeric(log2.counts[GeneNames[x], ])
  vals=vals[1:8]
  myDF=as.tibble(cbind(as.character(condition), as.numeric(vals)))
  colnames(myDF)=c("Condition", "Value")
  name=paste(FriendlyName[x], "Normalized Transcript Counts - Plasma 1 hr vs Control", "Log 2 Fold Change =",
             round(Genes[GeneNames[x],]$log2FoldChange, digits=2), "FDR Adjusted P Value =", scientific(Genes[GeneNames[x], ]$padj, digits=2))
  var=paste("plot", x, sep="")
  plot=ggplot(myDF, aes(x=as.numeric(Value), fill=as.character(Condition))) +
    geom_density(alpha=0.4) +
    xlab("Transcript Count") +
    ylab("Density") +
    ggtitle(name) + 
    labs(fill= "Condition") +
    theme(text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  #plotList=append(plotList, plot, after=length(plotList))
  assign(var, plot)
  #name=paste("densityPlot", FriendlyName[x], sep="-")
  #ggsave(name, plot = last_plot(), device = "png", dpi = 500, limitsize = FALSE)
}

library(cowplot)
#combined = plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, plot13, plot14, ncol=4, nrow=4)
combined=plot_grid(plot1, plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10, plot11, plot12, nrow=6, ncol=2)
ggsave("combined.png", plot=combined, device = "png", width = 25, height=45, dpi = 300, limitsize = FALSE)