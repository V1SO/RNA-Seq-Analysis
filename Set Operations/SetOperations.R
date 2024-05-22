library(venneuler)
one24intersect=intersect(rownames(res.sig.1hr), rownames(res.sig.24hr))
one24SetSize=(length(res.sig.1hr$baseMean)+length(res.sig.24hr$baseMean))
one24OverlapPercent=100*(length(one24intersect)/(one24SetSize/2))

one48intersect=intersect(rownames(res.sig.1hr), rownames(res.sig.48hrTMZ))
one48SetSize=(length(res.sig.1hr$baseMean)+length(res.sig.48hrTMZ$baseMean))
one48OverlapPercent=100*(length(one48intersect)/(one48SetSize/2))

one48TMZ24intersect=intersect(rownames(res.sig.24hr), rownames(res.sig.TMZ24))
one48TMZ24setDiff=setdiff(rownames(res.sig.TMZ24), rownames(res.sig.24hr))
one48TMZ24SetSize=(length(res.sig.24hr$baseMean)+length(res.sig.TMZ24$baseMean))
one48TMZ24overlapPercent=100*(length(one48TMZ24intersect)/(one48TMZ24SetSize/2))

TMZGenesSetDiff=res.sig.TMZ24[one48TMZ24setDiff,]
TMZGenesIntersection=res.sig.TMZ24[one48TMZ24intersect,]
write.csv(as.data.frame(TMZGenesSetDiff), file = "GenesThatAreDifferentiallyExpressedInTMZButNotInPlasma.csv")
write.csv(as.data.frame(TMZGenesIntersection), file = "GenesThatAreTheSameInTMZAndPlasma.csv")
index=rownames(TMZGenesSetDiff)[vikasList]
index=which(rownames(TMZGenesSetDiff)=="ENSG00000007952" |
              rownames(TMZGenesSetDiff)=="ENSG00000074771" |
              rownames(TMZGenesSetDiff)=="ENSG00000086991" |
              rownames(TMZGenesSetDiff)=="ENSG00000165168" |
              rownames(TMZGenesSetDiff)=="ENSG00000285441" |
              rownames(TMZGenesSetDiff)=="ENSG00000143878" |
              rownames(TMZGenesSetDiff)=="ENSG00000142168" |
              rownames(TMZGenesSetDiff)=="ENSG00000112096" |
              rownames(TMZGenesSetDiff)=="ENSG00000109572" |
              rownames(TMZGenesSetDiff)=="ENSG00000141837" |
              rownames(TMZGenesSetDiff)=="ENSG00000151067" |
              rownames(TMZGenesSetDiff)=="ENSG00000285479" |
              rownames(TMZGenesSetDiff)=="ENSG00000157388" |
              rownames(TMZGenesSetDiff)=="ENSG00000198216" |
              rownames(TMZGenesSetDiff)=="ENSG00000006283" |
              rownames(TMZGenesSetDiff)=="ENSG00000196557" |
              rownames(TMZGenesSetDiff)=="ENSG00000007402" |
              rownames(TMZGenesSetDiff)=="ENSG00000165995" |
              rownames(TMZGenesSetDiff)=="ENSG00000108878" |
              rownames(TMZGenesSetDiff)=="ENSG00000178363" |
              rownames(TMZGenesSetDiff)=="ENSG00000178372" |
              rownames(TMZGenesSetDiff)=="ENSG00000169885" |
              rownames(TMZGenesSetDiff)=="ENSG00000058404" |
              rownames(TMZGenesSetDiff)=="ENSG00000240583" |
              rownames(TMZGenesSetDiff)=="ENSG00000167580" |
              rownames(TMZGenesSetDiff)=="ENSG00000165272" |
              rownames(TMZGenesSetDiff)=="ENSG00000171885" |
              rownames(TMZGenesSetDiff)=="ENSG00000161798" |
              rownames(TMZGenesSetDiff)=="ENSG00000086159" |
              rownames(TMZGenesSetDiff)=="ENSG00000165269" |
              rownames(TMZGenesSetDiff)=="ENSG00000103375" |
              rownames(TMZGenesSetDiff)=="ENSG00000103569" |
              rownames(TMZGenesSetDiff)=="ENSG00000143595" |
              rownames(TMZGenesSetDiff)=="ENSG00000178301" |
              rownames(TMZGenesSetDiff)=="ENSG00000047936" |
              rownames(TMZGenesSetDiff)=="ENSG00000136888" |
              rownames(TMZGenesSetDiff)=="ENSG00000005381" |
              rownames(TMZGenesSetDiff)=="ENSG00000122986" |
              rownames(TMZGenesSetDiff)=="ENSG00000128340" |
              rownames(TMZGenesSetDiff)=="ENSG00000146648" |
              rownames(TMZGenesSetDiff)=="ENSG00000135679" |
              rownames(TMZGenesSetDiff)=="ENSG00000198625" |
              rownames(TMZGenesSetDiff)=="ENSG00000125995" |
              rownames(TMZGenesSetDiff)=="ENSG00000198015" |
              rownames(TMZGenesSetDiff)=="ENSG00000114026" |
              rownames(TMZGenesSetDiff)=="ENSG00000120129" |
              rownames(TMZGenesSetDiff)=="ENSG00000100292" |
              rownames(TMZGenesSetDiff)=="ENSG00000111206" )
candidateDiff=TMZGenesSetDiff[index,]
write.csv(as.data.frame(candidateDiff), file = "CandidateGenesThatAreDifferentiallyExpressedInTMZButNotInPlasma.csv")

index=which(rownames(TMZGenesIntersection)=="ENSG00000007952" |
              rownames(TMZGenesIntersection)=="ENSG00000074771" |
              rownames(TMZGenesIntersection)=="ENSG00000086991" |
              rownames(TMZGenesIntersection)=="ENSG00000165168" |
              rownames(TMZGenesIntersection)=="ENSG00000285441" |
              rownames(TMZGenesIntersection)=="ENSG00000143878" |
              rownames(TMZGenesIntersection)=="ENSG00000142168" |
              rownames(TMZGenesIntersection)=="ENSG00000112096" |
              rownames(TMZGenesIntersection)=="ENSG00000109572" |
              rownames(TMZGenesIntersection)=="ENSG00000141837" |
              rownames(TMZGenesIntersection)=="ENSG00000151067" |
              rownames(TMZGenesIntersection)=="ENSG00000285479" |
              rownames(TMZGenesIntersection)=="ENSG00000157388" |
              rownames(TMZGenesIntersection)=="ENSG00000198216" |
              rownames(TMZGenesIntersection)=="ENSG00000006283" |
              rownames(TMZGenesIntersection)=="ENSG00000196557" |
              rownames(TMZGenesIntersection)=="ENSG00000007402" |
              rownames(TMZGenesIntersection)=="ENSG00000165995" |
              rownames(TMZGenesIntersection)=="ENSG00000108878" |
              rownames(TMZGenesIntersection)=="ENSG00000178363" |
              rownames(TMZGenesIntersection)=="ENSG00000178372" |
              rownames(TMZGenesIntersection)=="ENSG00000169885" |
              rownames(TMZGenesIntersection)=="ENSG00000058404" |
              rownames(TMZGenesIntersection)=="ENSG00000240583" |
              rownames(TMZGenesIntersection)=="ENSG00000167580" |
              rownames(TMZGenesIntersection)=="ENSG00000165272" |
              rownames(TMZGenesIntersection)=="ENSG00000171885" |
              rownames(TMZGenesIntersection)=="ENSG00000161798" |
              rownames(TMZGenesIntersection)=="ENSG00000086159" |
              rownames(TMZGenesIntersection)=="ENSG00000165269" |
              rownames(TMZGenesIntersection)=="ENSG00000103375" |
              rownames(TMZGenesIntersection)=="ENSG00000103569" |
              rownames(TMZGenesIntersection)=="ENSG00000143595" |
              rownames(TMZGenesIntersection)=="ENSG00000178301" |
              rownames(TMZGenesIntersection)=="ENSG00000047936" |
              rownames(TMZGenesIntersection)=="ENSG00000136888" |
              rownames(TMZGenesIntersection)=="ENSG00000005381" |
              rownames(TMZGenesIntersection)=="ENSG00000122986" |
              rownames(TMZGenesIntersection)=="ENSG00000128340" |
              rownames(TMZGenesIntersection)=="ENSG00000146648" |
              rownames(TMZGenesIntersection)=="ENSG00000135679" |
              rownames(TMZGenesIntersection)=="ENSG00000198625" |
              rownames(TMZGenesIntersection)=="ENSG00000125995" |
              rownames(TMZGenesIntersection)=="ENSG00000198015" |
              rownames(TMZGenesIntersection)=="ENSG00000114026" |
              rownames(TMZGenesIntersection)=="ENSG00000120129" |
              rownames(TMZGenesIntersection)=="ENSG00000100292" |
              rownames(TMZGenesIntersection)=="ENSG00000111206" )
candidateIntersection=TMZGenesIntersection[index,]
write.csv(as.data.frame(candidateIntersection), file = "CandidateGenesThatAreTheSameInTMZAndPlasma.csv")

allIntersect=intersect(one24intersect,one48intersect)

commonGenes=cbind.data.frame(as.character(res.sig.1hr[allIntersect, ]$symbol),
                  as.numeric(res.sig.1hr[allIntersect, ]$padj),
                  as.numeric(res.sig.24hr[allIntersect, ]$padj),
                  as.numeric(res.sig.48hrTMZ[allIntersect, ]$padj),
                  as.numeric(res.sig.1hr[allIntersect, ]$log2FoldChange),
                  as.numeric(res.sig.24hr[allIntersect, ]$log2FoldChange),
                  as.numeric(res.sig.48hrTMZ[allIntersect, ]$log2FoldChange))
colnames(commonGenes)=c("Gene", "1 hr vs Control Adjusted P Value",
                        "24 hr vs Control Adjusted P Value",
                        "48 hr + TMZ vs Control Adjusted P Value", 
                        "1 hr vs Control Log Fold Change",
                        "24 hr vs Control Log Fold Change",
                        "48 hr + TMZ vs Control Log Fold Change")
commonGenes124=cbind.data.frame(as.character(res.sig.1hr[one24intersect, ]$symbol),
                             as.numeric(res.sig.1hr[one24intersect, ]$padj),
                             as.numeric(res.sig.24hr[one24intersect, ]$padj),
                             as.numeric(res.sig.1hr[one24intersect, ]$log2FoldChange),
                             as.numeric(res.sig.24hr[one24intersect, ]$log2FoldChange))
colnames(commonGenes124)=c("Gene", "1 hr vs Control Adjusted P Value",
                        "24 hr vs Control Adjusted P Value",
                        "1 hr vs Control Log Fold Change",
                        "24 hr vs Control Log Fold Change")
#commonGenes=data.frame(commonGenes)
variance=rowVars(commonGenes124[,c(4,5)])
commonGenes124$variance=variance
write.csv(as.data.frame(diffExpressed), file = "CandidateGenesDifferentiallyExpressed1hrand24hr.csv")