#Note: importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(biomaRt)
library(Biobase)
library(limma)
library(edgeR)
library(HTSFilter)
library(ggplot2)
library(BiocParallel)
library(gplots)
library(vsn)
library(ggrepel)

### Gather all counts files
setwd("~/tyson")
counts.data <- read.delim("countsallVikas.txt",sep = "\t",header = T,row.names = 1)
metastuff <- read.delim("metaData.txt",sep = "\t", header = T)
setwd("~/tyson/VikasDifferentialExpression24hrPlasmavs48hrsPlasma+TMZ")
colnames(counts.data) <- gsub("._Aligned_sorted","",colnames(counts.data))
colnames(counts.data) <- gsub("_ForwardPaired_","",colnames(counts.data))
colnames(counts.data) <- gsub(".._.","",colnames(counts.data))
colnames(counts.data) <- gsub(".fastq","",colnames(counts.data))

counts.data <- counts.data[,c(6:21)]
counts.data <- counts.data[,order(colnames(counts.data))]
#rownames(counts.data) <- substr(rownames(counts.data),1,15)

### Cross reference the counts data with sample meta data
counts.data <- counts.data[rowSums(counts.data) > 1, ] #Filter out those features that are only zeros
dim(counts.data) #30009 13
#metastuff <- metastuff[,order(metastuff$ID)]
colnames(counts.data) <- metastuff$X.Sample_title

###Build a DESeq object
#colnames(counts.data) <- c(1:ncol(counts.data)) #Unique names are required for creating DESeq objects, but will be substituted after technical replicate collapsing
colnames(metastuff)[colnames(metastuff)=="X.Sample_characteristics_ch1"] <- "Treatment"
metastuff[is.na(metastuff)] <- 0
metastuff$Treatment <- as.factor(metastuff$Treatment) #44- SLE, 44-CTL
#metastuff$WHO <- as.factor(metastuff$WHO)
dds <- DESeqDataSetFromMatrix(countData = counts.data, colData = metastuff, design = ~ Treatment)

gene.list <- as.data.frame(rownames(counts.data))
colnames(gene.list) <- "ensembl"
gene.list <- substr(gene.list$ensembl,1,15)
gene.list <- as.data.frame(gene.list)
colnames(gene.list) <- c("ensembl")
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
#                    host = "uswest.ensembl.org", 
#                    ensemblRedirect = FALSE )
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="uswest.ensembl.org", ensemblRedirect = FALSE)
genemap <- getBM(attributes = c("ensembl_gene_id_version","ensembl_gene_id","entrezgene_id","hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = gene.list$ensembl,
                 mart = ensembl)
idx <- match(gene.list$ensembl, genemap$ensembl_gene_id)
gene.list$entrez <- genemap$entrezgene[idx]
gene.list$symbol <- genemap$hgnc_symbol[idx]
gene.list$ensembl_just_id <- genemap$ensembl_gene_id[idx]
rownames(dds) <- substr(rownames(dds),1,15)
annotation <- gene.list[match(rownames(dds), gene.list$ensembl),]
all(rownames(dds) == annotation$ensembl)


### Filter raw Glomeruli reads with HTSFilter to eliminate uninformative and repetitive transcripts by Jaccard Index
idx <- as.character(dds$Treatment) #Make a vector of WHO variables for HTS filter
#setwd("~/Desktop/Roche_dataset")
pdf("HTS-Filter_Plot.pdf", width = 16, height = 9)
filtered <- HTSFilter(dds,
                      conds = idx,
                      pAdjustMethod = "BH",
                      s.len = 100,
                      s.min = 1,
                      s.max = 200,
                      normalization = "DESeq",
                      plot = T)
dev.off()
filtered <- filtered$filteredData #Get the data frame of the filtered transcript reads
dim(filtered) #11756 29 remaining features
filtered@rowRanges@elementMetadata@listData$ensembl <- droplevels(filtered@rowRanges@elementMetadata@listData$ensembl)

#dds_SLE_CTL <- dds[,dds$Treatment=="CTL" | dds$Treatment=="SLE"] 
#dds_pSS_CTL <- dds[,dds$Treatment=="HD" | dds$Treatment=="pSS"]


dds <- DESeq(filtered) #Normalizes via DESeq2 method; already done with HTSFilter, but results will be the same
pdf(file="Dispersion_Estimates.pdf",height = 9,width = 16)
plotDispEsts( dds, ylim = c(1e-8, 1e1) ) #Check dispersion estimates 
dev.off()


### Differential Analyses
## Compare the SLE versus CTL
#dds <- dds[,dds$Treatment=="CTL" | dds$Treatment=="SLE"] 
res <- results(dds,
               contrast=c("Treatment", "Plasma48hrTMZ", "Plasma24hr"),
               alpha = 0.05,
               pAdjustMethod = "BH",
               independentFiltering=FALSE,  #Since HTSFilter was used, no need for DESeq to filter
               BPPARAM = SnowParam(4),
               parallel = T)
res <- res[order(res$padj),]
gene.ids <- annotation[match(rownames(res), annotation$ensembl),]
res$symbol <- gene.ids$symbol
res$entrez <- gene.ids$entrez
index <- which(res$padj<0.05)
res.sig <- res[index,]
res.sig <- as.data.frame(res.sig)
res.sig <- res.sig[ order(res.sig$padj), ] 
res.sig <- res.sig[!duplicated(res.sig$entrez), ]
res.sig <- res.sig[!is.na(res.sig$entrez), ]
#setwd("~/Desktop/Roche_dataset")
write.csv(as.data.frame(res.sig), file = "RNAseq_DESeq_SLE-CTL_TSD_29October2019.csv")
dim(res.sig) #1197 DE genes


### Check MA plots of differential data
#setwd("~/Dropbox (AMPEL BioSolutions)/Bioinformatics/studies_sle/adult/GSE112087_RNAseq_wholeblood/Deseq2")
pdf(file="MA_plot.pdf",height = 9,width = 16)
DESeq2::plotMA( dds, ylim = c(-10, 10), main = "All genes")
dev.off()


### Data Transformations
## Glomeruli
# Log2 transformation
log2 <- normTransform(dds)
log2.counts <- assay(log2)
log2.counts <- as.data.frame(log2.counts)
## Alternatively:
# glom.log2.counts <- log2(1+counts(ddssva.glom,normalized=TRUE))
# Regularized log transformation
rld <- rlogTransformation(dds, blind=FALSE)
rld.counts <- assay(rld)
rld.counts <- as.data.frame(rld.counts)
# Variance stablizing transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.counts <- assay(vsd)
vsd.counts <- as.data.frame(vsd.counts)

gene.ids <- annotation[match(rownames(log2.counts), annotation$ensembl),]
log2.counts$symbol <- gene.ids$symbol
log2.counts$entrez <- gene.ids$entrez
#log2.counts <- log2.counts[ ,c((ncol(log2.counts)-1):ncol(log2.counts),1:(ncol(log2.counts)-2))]
write.csv(as.data.frame(log2.counts), file = "DESeq_log2-counts_TSD_Vikas_December_2019.csv")
rld.counts$symbol <- gene.ids$symbol
rld.counts$entrez <- gene.ids$entrez
#rld.counts <- rld.counts[ ,c((ncol(glom.rld.counts)-1):ncol(glom.rld.counts),1:(ncol(glom.rld.counts)-2))]
write.csv(as.data.frame(rld.counts), file = "DESeq_rlog-counts_TSD_Vikas_December_2019.csv")
vsd.counts$symbol <- gene.ids$symbol
vsd.counts$entrez <- gene.ids$entrez
#glom.vsd.counts <- glom.vsd.counts[ ,c((ncol(glom.vsd.counts)-1):ncol(glom.vsd.counts),1:(ncol(glom.vsd.counts)-2))]
write.csv(as.data.frame(vsd.counts), file = "DESeq_VST-counts_TSD_Vikas_December_2019.csv")


### Compare Data Transformation Methods
pdf(file="DT_Transformation_Comparison_1.pdf",height=9,width=16)
par(mfrow=c(1,3))
plot(assay(log2)[,1:4],col="#00000020",pch=20,cex=0.3,main="Log2 Transformed")
plot(assay(rld)[,1:4],col="#00000020",pch=20,cex=0.3,main="Regularized Logarithm (rLog)")
plot(assay(vsd)[,1:4],col="#00000020",pch=20,cex=0.3,main="Variance Stabilization Tranformation (VST)")
dev.off()

pdf(file="DT_Transformation_Comparison_2.pdf",height = 9,width = 16)
select <- order(rowMeans(assay(log2)),decreasing=TRUE)[1:30]
colors <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(assay(log2)[select,], col = colors,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(9,16), main = "Log2 Transformed")
heatmap.2(assay(rld)[select,], col = colors,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(9,16), main = "Regularized Logarithm (rLog)")
heatmap.2(assay(vsd)[select,], col = colors,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(9,16), main = "Variance Stabilization Tranformation (VST)")
dev.off()

pdf(file="DT_Transformation_Comparison_3.pdf",height=9,width=16)
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(assay(log2[notAllZero,]))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()

# Show flattening of data toward consistent variance
pdf(file="DT_Variance_Flattening.pdf",height = 9,width = 16)
par(mai=ifelse(1:4 <= 2, par('mai'), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord]<150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c('blue', 'red')
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type="p", lty=1, col=vstcol, xlab='n', ylab='f(n)')
legend('bottomright', legend = c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

# Consider the total distance between samples
pdf("DT_PCA_rLog.pdf",height = 9, width = 16)
print( plotPCA( rld, intgroup = c( "Treatment")))
dev.off()

pdf(file="DT_Sample_Distances_rLog.pdf",height = 9,width = 16)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld@colData@rownames
colnames(sampleDistMatrix) <- rld@colData@rownames
hc <- hclust(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
heatmap.2(sampleDistMatrix,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=colors,margins=c(9,16),
          main = "Regularized Logarithm (rLog)")
dev.off()

pdf("DT_PCA_VST.pdf",height = 9, width = 16)
print( plotPCA( vsd, intgroup = c( "Treatment")) )
dev.off()

pdf(file="DT_Sample_Distances_VST.pdf",height = 9,width = 16)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd@colData@rownames
colnames(sampleDistMatrix) <- vsd@colData@rownames
hc <- hclust(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
heatmap.2(sampleDistMatrix,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=colors,margins=c(9,16),
          main = "Variance Stabilization Tranformation (VST)")
dev.off()

pdf("DT_PCA_Log2.pdf",height = 9, width = 16)
print( plotPCA( log2, intgroup = c("Treatment")))
dev.off()

pdf(file="DT_Sample_Distances_Log2.pdf",height = 9,width = 16)
sampleDists <- dist(t(assay(log2)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- log2@colData@rownames
colnames(sampleDistMatrix) <- log2@colData@rownames
hc <- hclust(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
heatmap.2(sampleDistMatrix,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=colors,margins=c(9,16),
          main = "Log2")
dev.off()

# Display the DEGs
#topVarGenes <- order( rowVars( assay(log2[which(rownames(log2)%in%rownames(res.sig))]) ), decreasing=TRUE )#, 100 )
#pdf("Most_Variant_Genes.pdf",height = 9,width = 16)
#heatmap.2( assay(log2)[ topVarGenes, ], scale="row",
#           trace="none", dendrogram="both",
#           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
#           margins = c(16,9))
#dev.off()


###Save results
save(dds,
     dds,
     res,
     res.sig,
     log2,
     log2.counts,
     rld,
     rld.counts,
     vsd,
     vsd.counts,
     file="VikasDEseq.RData")

##### Make an ESET for Glomeruli
# Gather and organize pData, fData, and expression data
pData <- metastuff
rownames(pData) <- pData$X.Sample_title
#pData <- pData[order(rownames(pData)), ]
exprs <- log2.counts[,c(1:16)]
exprs <- exprs[,order(colnames(exprs))]
identical(colnames(exprs),rownames(pData))
fData <- as.data.frame(log2.counts[,c(17:18)])
fData$ensembl <- rownames(fData)
exprs <- exprs[order(rownames(exprs)), ]
fData <- fData[(order(fData$ensembl)),]
identical(rownames(fData),rownames(exprs))

###########################################################
# CREATE ESET
# Biobase expects a matrix of expression values to generate an eset
exprs.eset <- as.matrix(exprs)
rownames(exprs.eset) <- rownames(exprs)

pData.meta <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData) )
pData.eset <- new("AnnotatedDataFrame", data = data.frame(pData), varMetadata = pData.meta)

fData.meta <- data.frame(labelDescription = colnames(fData), row.names = colnames(fData) )
fData.eset <- new("AnnotatedDataFrame", data = data.frame(fData), varMetadata = fData.meta)

expData <- new("MIAME", name = "PB",
               lab = "Roche_GT", 
               contact = "PB",title = "Human CD4T RNA-seq",
               abstract = "blah,blah,blah",
               url = "ampel.com",
               other = list(notes = "stuff and things"));

eset <- ExpressionSet(assayData=exprs.eset,
                      phenoData = pData.eset,
                      featureData = fData.eset,
                      experimentData = expData,
                      annotation = "Human CD4T RNA-seq HTS Filtered 26April2019")

eset
table(eset$Treatment)
save(eset, file="VikasTysonEset.RData")
# assign up or down
res.sig$fc_direction <- "UP"
index <- which(res.sig$log2FoldChange<=0)
res.sig[ index, "fc_direction"] <- "DOWN"
#Create scatterplot
reallysig=as.data.frame(res.sig[1:5000,])
res.sig.TMZ24=res.sig

pdf(file="DifferentiallyExpressed.pdf",height = 10,width = 10)

ggplot(reallysig, aes(x=-log(padj), y=log2FoldChange)) + geom_point() +
  geom_label(label=reallysig$symbol, nudge_y = 0.1)

dev.off()


pdf(file="LabelledDifferentiallyExpressed2.pdf",height = 12,width = 12)
genes <- reallysig
genes$Significant <- ifelse(genes$padj < 4*(10^-300), "FDR < 4e-300", "Not Sig")
ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant, size=baseMean)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 20) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(genes, padj < 4*(10^-300)),
    aes(label = symbol),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
dev.off()

pdf(file="zscore.pdf",height = 40,width = 40)
genes <- res.sig.TMZ24
genes$Significant <- ifelse(genes$padj < 4*(10^-250), "FDR < 4e-250", "Not Sig")
plot=ggplot(genes, aes(x = log2FoldChange, y = stat)) +
  geom_point(aes(color = Significant, size=baseMean)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 20) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(genes, padj < 4*(10^-250)),
    aes(label = symbol),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
plot+ylim(-400, 400)
dev.off()
