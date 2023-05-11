rm(list=ls())
library(org.Mm.eg.db)
library(edgeR)
library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(pheatmap)
library(AnnotationDbi)  
library(ChIPpeakAnno)
library(tweeDEseq)
library(tidyverse)
library(MACSr)
library(Herper)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DiffBind)


# RNA-Seq

rawCTCMP1 <- read.table("~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/CMP/ENCFF623OLU.tsv", header = TRUE, row.names = 1)
rawCTCMP2 <- read.table("~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/CMP/ENCFF691MHW.tsv", header = TRUE, row.names = 1)
rawCTEry1 <- read.table("~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/Erythroblast/ENCFF858JHF.tsv", header = TRUE, row.names = 1)
rawCTEry2 <- read.table("~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/Erythroblast/ENCFF342WUL.tsv", header = TRUE, row.names = 1)
rawCTHSC1 <- read.table("~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/HSC/ENCFF064MKY.tsv", header = TRUE, row.names = 1)
rawCTHSC2 <- read.table("~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/HSC/ENCFF247FEJ.tsv", header = TRUE, row.names = 1)



## RNA-seq CMP(2) vs Erythroblast(4)

### DESeq2

####Data Prep and Cleaning
id<- rawCTEry1$transcript_id.s.
CMP1ExpCounts <- rawCTCMP1$expected_count
CMP2ExpCounts <- rawCTCMP2$expected_count
Ery1ExpCounts <- rawCTEry1$expected_count
Ery2ExpCounts <- rawCTEry2$expected_count
Counts24 <- cbind(id,CMP1ExpCounts,CMP2ExpCounts,Ery1ExpCounts,Ery2ExpCounts)
rownames(Counts24) <- id
Counts24F<-Counts24[,-1]
Counts24F<-data.matrix(Counts24F)
class(Counts24F)<-"numeric"
colnames(Counts24F)
colData<-data.frame("CellType" = c("CMP","CMP","Ery","Ery"))
rownames(colData)<-colnames(Counts24F)




####Data set into DE-Seq
dds0 <- DESeqDataSetFromMatrix(countData = round(Counts24F), colData = colData, design= ~ CellType)
keep <- rowSums(counts(dds0)) >= 10 # keep transcripts with more than 10 reads
dds <- dds0[keep,]
dds$CellType <- factor(dds$CellType, levels =c("CMP","Ery"))
dds<-cpm(dds)


####Des-seq results
dds.result <- DESeq(dds)
res <- results(dds.result)
resOrdered <- res[order(res$pvalue),] # order results by the smallest p-value
DESQ <- as.data.frame(res)
summary(res)
sum(res$padj < 0.01, na.rm=TRUE)



####Dispersion Analysis
dds1 <- estimateSizeFactors(dds) 
dds2 <- estimateDispersions(dds1) 
plotDispEsts(dds2, cex.lab=1.5 ) 
dds3 <- nbinomWaldTest(dds2) 
res2 <- results(dds3)
hist(res2$pvalue, xlab='Raw p-value', cex.lab=1.5, cex.main=2)

####significant genes
t <- res$padj < .01 & !is.na(res$padj)
sig.gn <- rownames(assay(dds))[t]
length(sig.gn)



### Limma-voom
#### prepare data
d0 <- DGEList(round(Counts24F),group=rep(c("HSC","EBT"),each=2))
d0 <- calcNormFactors(d0)
cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)



####create model matrix
group <- c("CMP","CMP","Ery","Ery")
mm <- model.matrix(~group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)


####fit matrix model
fit <- lmFit(y, mm)
head(coef(fit))


###make contrasts
contr <- makeContrasts(-groupEry, levels = colnames(coef(fit)))
row.names(contr)<-c("(Intercept)","groupEry")
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.01))

### plot MA plot from limma-voom
limma::plotMA(fit, coef=2, main="CMP vs Ery comparison")
abline(h=0, col="red", lwd=2)
hist(top.table$adj.P.Val,xlab="Adjusted Pvalues",main="Hist of P.adj")


### Enrichment
library(biomaRt)
library(org.Mm.eg.db)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
transid<-sig.gn
resont <- getBM(attributes = c('ensembl_transcript_id_version', 
                               'ensembl_gene_id', 
                               'external_transcript_name',
                               'external_gene_name'),
                filters = 'ensembl_transcript_id_version', 
                values = transid,
                mart = mart)
GORNA1<-clusterProfiler::enrichGO(resont$ensembl_gene_id,keyType = "ENSEMBL" ,OrgDb = 'org.Mm.eg.db', ont="BP")
goplot(GORNA1)
barplot(GORNA1, showCategory=20)
enrichplot::upsetplot(GORNA1)


### Clustering & Heatmaps

library("pheatmap")
ntd <- normTransform(dds.result)
select <- order(rowMeans(counts(dds.result,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds.result)$CellType)
colnames(df) <- "Cell Type"
rownames(df) <- c("ENCFF623OLU", "ENCFF691MHW", "ENCFF858JHF", "ENCFF342WUL")
labs <- assay(ntd)[select,]
colnames(labs)<-c("ENCFF623OLU", "ENCFF691MHW", "ENCFF858JHF", "ENCFF342WUL")


pheatmap(labs,
         show_rownames=FALSE,
         cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
         cluster_cols = TRUE,
         annotation_col=df)




# CMP(1) vs HSC(2)


### DESeq2

id<- rawCTEry1$transcript_id.s.
CMP1ExpCounts <- rawCTCMP1$expected_count
CMP2ExpCounts <- rawCTCMP2$expected_count
HSC1ExpCounts <- rawCTHSC1$expected_count
HSC2ExpCounts <- rawCTHSC2$expected_count
Counts24.2 <- cbind(id,CMP1ExpCounts,CMP2ExpCounts,HSC1ExpCounts,HSC2ExpCounts)
rownames(Counts24.2)<-id
Counts24F.2<-Counts24.2[,-1]
Counts24F.2<-data.matrix(Counts24F)
class(Counts24F.2)<-"numeric"
colnames(Counts24F.2)
colData.2<-data.frame("CellType" = c("CMP","CMP","HSC","HSC"))
rownames(colData.2)<-colnames(Counts24F.2)




dds0.2 <- DESeqDataSetFromMatrix(countData = round(Counts24F.2), colData = colData.2, design= ~ CellType)
keep.2 <- rowSums(counts(dds0.2)) >= 10 # keep transcripts with more than 10 reads
dds.2 <- dds0.2[keep.2,]
dds.2$CellType <- factor(dds.2$CellType, levels =c("CMP","HSC"))
dds.2<-cpm(dds.2)



dds.result.2 <- DESeq(dds.2)
res.2 <- results(dds.result.2)
resOrdered.2 <- res.2[order(res.2$pvalue),] # order results by the smallest p-value
DESQ.2 <- as.data.frame(res.2)
sum(res.2$padj < 0.01, na.rm=TRUE)




dds1.2 <- estimateSizeFactors(dds.2) 
dds2.2 <- estimateDispersions(dds1.2) 
plotDispEsts(dds2.2, cex.lab=1.5 ) 
dds3.2 <- nbinomWaldTest(dds2.2) 
res2.2 <- results(dds3.2)
hist(res2.2$pvalue, xlab='Raw p-value', cex.lab=1.5, cex.main=2)
t.2 <- res.2$padj < .01 & !is.na(res.2$padj)
sig.gn.2 <- rownames(assay(dds.2))[t.2]
length(sig.gn.2)
plotMA(dds.result.2)

### Enrichment

library(biomaRt)
library(org.Mm.eg.db)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
transid2<-sig.gn.2
resont2 <- getBM(attributes = c('ensembl_transcript_id_version', 
                                'ensembl_gene_id', 
                                'external_transcript_name',
                                'external_gene_name'),
                 filters = 'ensembl_transcript_id_version', 
                 values = transid2,
                 mart = mart)
GORNA2<-clusterProfiler::enrichGO(resont2$ensembl_gene_id,keyType = "ENSEMBL" ,OrgDb = 'org.Mm.eg.db', ont="BP")
goplot(GORNA2)
barplot(GORNA2, showCategory=20)
enrichplot::upsetplot(GORNA2)



### Clustering

library("pheatmap")
ntd.2 <- normTransform(dds.result.2)
select.2 <- order(rowMeans(counts(dds.result.2,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
df.2 <- as.data.frame(colData(dds.result.2)$CellType)
colnames(df.2) <- "Cell Type"
rownames(df.2) <- c("ENCFF623OLU", "ENCFF691MHW", "ENCFF064MKY", "ENCFF247FEJ")
labs.2 <- assay(ntd.2)[select.2,]
colnames(labs.2)<-c("ENCFF623OLU", "ENCFF691MHW", "ENCFF064MKY", "ENCFF247FEJ")


pheatmap(labs.2,
         show_rownames=FALSE,
         cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
         cluster_cols = TRUE,
         annotation_col=df.2)



# ATAC-Seq

library(ChIPseeker)
CMPvsEry <- dba(sampleSheet="~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/samplesheet1.csv")
CMPvsHSC <- dba(sampleSheet="~/Desktop/PSU/SeniorYear/SpringSemester/STAT555/FinalProjData/Reads/DATAwithBedFiles/samplesheet2.csv")
CMPvsEry.BL <- dba.blacklist(CMPvsEry,blacklist=DBA_BLACKLIST_MM10)
CMPvsHSC.BL <- dba.blacklist(CMPvsHSC,blacklist=DBA_BLACKLIST_MM10)

## Create counts

CMPvsEryCounts <- dba.count(CMPvsEry, summits=250)
CMPvsHSCCounts <- dba.count(CMPvsHSC, summits=250)


## CMP vs Ery
CMPvsEryCont <- dba.contrast(CMPvsEryCounts, categories=DBA_TREATMENT, minMembers = 2)
DESeqCMPvsEry <- dba.analyze(CMPvsEryCont, method=DBA_DESEQ2)
ResDESeqCMPvsEry <- dba.report(DESeqCMPvsEry,method=DBA_DESEQ2)

dba.plotMA(DESeqCMPvsEry,method=DBA_DESEQ2)
dba.plotHeatmap(DESeqCMPvsEry)

out <- as.data.frame(ResDESeqCMPvsEry)
write.table(out, file="CMPvsEry_deseq2.txt", sep="\t", quote=F, row.names=F)

CMPvsEry_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)
nrow(CMPvsEry_enrich)

# Write to file
write.table(CMPvsEry_enrich, file="CMPvsEry_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

## CMP vs HSC

CMPvsHSCCont <- dba.contrast(CMPvsHSCCounts, categories=DBA_TREATMENT, minMembers = 2)
DESeqCMPvsHSC <- dba.analyze(CMPvsHSCCont, method=DBA_DESEQ2)
v <- dba.report(DESeqCMPvsHSC,method=DBA_DESEQ2)

dba.plotMA(DESeqCMPvsHSC,method=DBA_DESEQ2)
dba.plotHeatmap(DESeqCMPvsHSC)

out2 <- as.data.frame(ResDESeqCMPvsHSC)
write.table(out2, file="CMPvsHSC_deseq2.txt", sep="\t", quote=F, row.names=F)
CMPvsHSC_enrich <- out2 %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)
nrow(CMPvsHSC_enrich)

# Write to file
write.table(CMPvsHSC_enrich, file="CMPvsHSC_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
### Clustering

counts_CMPvsHSC<- dba.peakset(DESeqCMPvsHSC, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
CMPvsHSC_enrichHM <- out2 %>% arrange(FDR)%>%
  filter(FDR < 0.05 & Fold < 0)
colnames(CMPvsHSC_enrichHM)<-c("CHR", "START", "END")
RawReads_CMPvsHSC<-merge(counts_CMPvsHSC,CMPvsHSC_enrichHM,by=c("CHR", "START", "END"))
RawReads_CMPvsHSC<-RawReads_CMPvsHSC[,4:7]
HCLUSTMatCMPvsHSC <-dist(as.matrix(RawReads_CMPvsHSC[1:5000,]))
dendoMatCMPvsEry <- hclust(HCLUSTMatCMPvsHSC,method = "complete")
plot(dendoMatCMPvsEry, labels = FALSE)

####Heat maps
library("pheatmap")
df.3 <- DESeqCMPvsHSC$samples
df.3 <-df.3[,1:2]
rownames(df.3) <- df.3$SampleID
df.3$CellType<-df.3$Treatment
df.3$SampleID<-NULL
df.3$Treatment<-NULL
labs.3 <- as.matrix(RawReads_CMPvsHSC[1:50,])
colnames(labs.3)<-CMPvsHSCCountsClus$samples$SampleID


pheatmap(labs.3,
         show_rownames=FALSE,
         cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
         cluster_cols = TRUE,
         annotation_col=df.3)

####clustering
counts_CMPvsEry<- dba.peakset(DESeqCMPvsEry, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
CMPvsEry_enrichHM <- out  %>% arrange(FDR)%>%
  filter(FDR < 0.05 & Fold < 0)
colnames(CMPvsEry_enrichHM)<-c("CHR", "START", "END")
RawReads_CMPvsEry<-merge(counts_CMPvsEry,CMPvsEry_enrichHM,by=c("CHR", "START", "END"))
RawReads_CMPvsEry<-RawReads_CMPvsEry[,4:7]
HCLUSTMatCMPvsEry <-dist(as.matrix(RawReads_CMPvsEry[1:500,]))
dendoMatCMPvsEry <- hclust(HCLUSTMatCMPvsEry,method = "complete")
plot(dendoMatCMPvsEry, labels = FALSE)

###ATAC-seq Heat maps
library("pheatmap")
df.4 <- DESeqCMPvsEry$samples
df.4 <-df.4[,1:2]
rownames(df.4) <- df.4$SampleID
df.4$CellType<-df.4$Treatment
df.4$SampleID<-NULL
df.4$Treatment<-NULL
labs.4 <- as.matrix(na.omit(RawReads_CMPvsEry[1:50,]))
colnames(labs.4)<-DESeqCMPvsEry$samples$SampleID

pheatmap(labs.4,
         show_rownames=FALSE,
         cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
         cluster_cols = TRUE,
         annotation_col=df.4)


## Enrichment Analysis
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakCMPvsEry <- readPeakFile("CMPvsEry_enriched.bed")
peakCMPvsHSC <- readPeakFile("CMPvsHSC_enriched.bed")

peakAnnoCMPvsEry <- annotatePeak(peakCMPvsEry, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoCMPvsHSC <- annotatePeak(peakCMPvsHSC, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
###plots from ATAC-seq
plotAnnoPie(peakAnnoCMPvsEry)
plotAnnoPie(peakAnnoCMPvsHSC)

library(ggupset)
enrichplot::upsetplot(peakAnnoCMPvsEry)
enrichplot::upsetplot(peakAnnoCMPvsHSC)

###pathway enrichment
library(ReactomePA)
library(enrichplot)
pathway1 <- enrichPathway(as.data.frame(peakAnnoCMPvsEry)$geneId,organism = "mouse")
head(pathway1, 2)
gene <- seq2gene(as.GRanges(peakAnnoCMPvsEry),tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene,organism = "mouse")
head(pathway2, 2)
dotplot(pathway2)
###gene enrichment
library(enrichplot)
clusterPro_CMPvsEry <- clusterProfiler::enrichGO(as.data.frame(peakAnnoCMPvsEry)$geneId, OrgDb = 'org.Mm.eg.db', ont="BP")
goplot(clusterPro_CMPvsEry)
barplot(clusterPro_CMPvsEry, showCategory=20)
enrichplot::upsetplot(clusterPro_CMPvsEry)

###pathway enrichment
library(ReactomePA)
library(enrichplot)
pathway1 <- enrichPathway(as.data.frame(peakAnnoCMPvsHSC)$geneId,organism = "mouse")
head(pathway1, 2)
gene <- seq2gene(as.GRanges(peakAnnoCMPvsHSC),tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene,organism = "mouse")
head(pathway2, 2)
dotplot(pathway2)
###gene enrichment
library(enrichplot)
clusterPro_CMPvsHSC <- clusterProfiler::enrichGO(as.data.frame(peakAnnoCMPvsHSC)$geneId, OrgDb = 'org.Mm.eg.db', ont="BP")
goplot(clusterPro_CMPvsHSC)
barplot(clusterPro_CMPvsHSC, showCategory=20)
enrichplot::upsetplot(clusterPro_CMPvsHSC)


