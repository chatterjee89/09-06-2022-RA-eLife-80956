library(edgeR)
library(gprofiler2)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(reticulate)
library(ggplot2)
library(org.Dm.eg.db)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(gplots)
library(DESeq2)
library(lazyeval)
library(knitr)

targets<-read.csv("/home/dchatterjee/Documents/bulkRNASeq/XFW_bulkRNASeqtumor_results/sample_info_tumor.csv")
group<-factor(paste(targets$Genotype,targets$Sex,sep="."))
cbind(targets,group=group)
x<-read.csv("/home/dchatterjee/Documents/bulkRNASeq/XFW_bulkRNASeqtumor_results/tumor_gene_counts.csv", row.names="Geneid")
cat("Number of input rows:\n")
nrow(x)
genetable=data.frame(gene_id=rownames(x))
y<-DGEList(counts=x, group=group, genes=genetable)
head(x)
y$samples
cat ("Counts:\n")
head(y$counts)
cat ("Genes:\n")
head(y$genes)
cperm <- cpm(y)
cat("Summary of Counts per Million:\n")
Summary<-summary(cperm)
write.csv(Summary, file="/home/dchatterjee/Documents/CPMSummary.csv")
countcheck <- cperm > 1
cat("After filtering for CPM < cut-off:\n")
head(countcheck)
keep<-which(rowSums(countcheck) >= 2)
y<-y[keep, ]
cat("Summary:\n")
summary<-summary(cpm(y))
write.csv(summary, file="/home/dchatterjee/Documents/FilteringSummary.csv")
logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:125]
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
tiff("Documents/heatmap_top100.tiff", height = 7000, width = 5000, res = 800)
pheatmap(highly_variable_lcpm, col=colorRampPalette(rev(brewer.pal(n = 10, name ="Spectral")))(100),  main="Top 100 most variable genes across samples",kmeans_k=NA, fontsize=4, cluster_rows=TRUE, cluster_cols=TRUE,  clustering_distance_rows = "euclidean", cutree_cols=7)
dev.off()
cat("Compared to the original summary:\n")
dim(y)
y<-calcNormFactors (y,method='TMM')
cat("Creating MDS plot: plotmDS.tiff.\n")
plotMDS(y)
cat("Done!\n")
design<-model.matrix(~0+group, data=y$samples)
colnames(design)<-levels(y$samples$group)
design
y<-estimateGLMCommonDisp(y, design)
y<-estimateGLMTrendedDisp(y, design)
y<-estimateGLMTagwiseDisp(y, design)
cat("Creating BCV plot: plotBCV.tiff\n")
plotBCV(y)
cat("Done!\n")
fit<-glmQLFit (y, design)
lrt<-glmLRT(fit)
summary(decideTests(lrt))
cat("Creating MD plot: plotMD.tiff\n")
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

my.contrasts1<-makeContrasts(WT_Time=WT.24H-WT.96H,
							Lgl_Time=Lgl_RNAi.24H-Lgl_RNAi.96H,
							Upd_Time=UAS_Upd.24H-UAS_Upd.96H,
							NICD_Time=UAS_NICD.24H-UAS_NICD.96H,
							UpdLgl_Time=Lgl_RNAi_UAS_Upd.24H-Lgl_RNAi_UAS_Upd.96H,
							NICDLgl_Time=Lgl_RNAi_UAS_NICD.24H-Lgl_RNAi_UAS_NICD.96H,
							levels=design)
qlf.WT_Time<-glmQLFTest(fit, contrast=my.contrasts1[,"WT_Time"])
qlf.Lgl_Time<-glmQLFTest(fit, contrast=my.contrasts1[,"Lgl_Time"])
qlf.Upd_Time<-glmQLFTest(fit, contrast=my.contrasts1[,"Upd_Time"])
qlf.NICD_Time<-glmQLFTest(fit, contrast=my.contrasts1[,"NICD_Time"])
qlf.UpdLgl_Time<-glmQLFTest(fit, contrast=my.contrasts1[,"UpdLgl_Time"])
qlf.NICDLgl_Time<-glmQLFTest(fit, contrast=my.contrasts1[,"NICDLgl_Time"])
qlf.WT_Time<-glmTreat(fit, contrast=my.contrasts1[,"WT_Time"], lfc=log2(1.8))
qlf.Lgl_Time<-glmTreat(fit, contrast=my.contrasts1[,"Lgl_Time"], lfc=log2(1.8))
qlf.Upd_Time<-glmTreat(fit, contrast=my.contrasts1[,"Upd_Time"], lfc=log2(1.8))
qlf.NICD_Time<-glmTreat(fit, contrast=my.contrasts1[,"NICD_Time"], lfc=log2(1.8))
qlf.UpdLgl_Time<-glmTreat(fit, contrast=my.contrasts1[,"UpdLgl_Time"], lfc=log2(1.8))
qlf.NICDLgl_Time<-glmTreat(fit, contrast=my.contrasts1[,"NICDLgl_Time"], lfc=log2(1.8))
de1<-decideTestsDGE(qlf.WT_Time, adjust.method="BH", p.value=0.05)
sum1<-summary(de1)
write.csv(sum1, file="/home/dchatterjee/Documents/Summary_WT_Time.csv")
de1<-decideTestsDGE(qlf.Lgl_Time, adjust.method="BH", p.value=0.05)
sum1<-summary(de1)
write.csv(sum1, file="/home/dchatterjee/Documents/Summary_Lgl_Time.csv")
de1<-decideTestsDGE(qlf.Upd_Time, adjust.method="BH", p.value=0.05)
sum1<-summary(de1)
write.csv(sum1, file="/home/dchatterjee/Documents/Summary_Upd_Time.csv")
de1<-decideTestsDGE(qlf.NICD_Time, adjust.method="BH", p.value=0.05)
sum1<-summary(de1)
write.csv(sum1, file="/home/dchatterjee/Documents/Summary_NICD_Time.csv")
de1<-decideTestsDGE(qlf.UpdLgl_Time, adjust.method="BH", p.value=0.05)
sum1<-summary(de1)
write.csv(sum1, file="/home/dchatterjee/Documents/Summary_UpdLgl_Time.csv")
de1<-decideTestsDGE(qlf.NICDLgl_Time, adjust.method="BH", p.value=0.05)
sum1<-summary(de1)
write.csv(sum1, file="/home/dchatterjee/Documents/Summary_NICDLgl_Time.csv")
tab1 <- topTags(qlf.WT_Time, n=Inf, sort.by="PValue", p.value=0.05)
write.csv(tab1, file="/home/dchatterjee/Documents/DGEgenes_WT_Time.csv")
tab1 <- topTags(qlf.Lgl_Time, n=Inf, sort.by="PValue", p.value=0.05)
write.csv(tab1, file="/home/dchatterjee/Documents/DGEgenes_Lgl_Time.csv")
tab1 <- topTags(qlf.Upd_Time, n=Inf, sort.by="PValue", p.value=0.05)
write.csv(tab1, file="/home/dchatterjee/Documents/DGEgenes_Upd_Time.csv")
tab1 <- topTags(qlf.NICD_Time, n=Inf, sort.by="PValue", p.value=0.05)
write.csv(tab1, file="/home/dchatterjee/Documents/DGEgenes_NICD.csv")
tab1 <- topTags(qlf.UpdLgl_Time, n=Inf, sort.by="PValue", p.value=0.05)
write.csv(tab1, file="/home/dchatterjee/Documents/DGEgenes_UpdLgl_Time.csv")
tab1 <- topTags(qlf.NICDLgl_Time, n=Inf, sort.by="PValue", p.value=0.05)
write.csv(tab1, file="/home/dchatterjee/Documents/DGEgenes_NICDLgl_Time.csv")
tiff("Documents/WT_Time.tiff", height = 3000, width = 3000, res = 400)
plotMD(qlf.WT_Time, main="WT 24H vs 96H", cex=0.75) + abline(h=c(-1, 1), col="blue")
dev.off()
tiff("Documents/Lgl_Time.tiff", height = 3000, width = 3000, res = 400)
plotMD(qlf.Lgl_Time, main="LglIRDcr2 24H vs 96H", cex=0.75) + abline(h=c(-1, 1), col="blue")
dev.off()
tiff("Documents/Upd_Time.tiff", height = 3000, width = 3000, res = 400)
plotMD(qlf.Upd_Time, main="Upd 24H vs 96H", cex=0.75) + abline(h=c(-1, 1), col="blue")
dev.off()
tiff("Documents/NICD_Time.tiff", height = 3000, width = 3000, res = 400)
plotMD(qlf.NICD_Time, main="NICD 24H vs 96H", cex=0.75) + abline(h=c(-1, 1), col="blue")
dev.off()
tiff("Documents/UpdLgl_Time.tiff", height = 3000, width = 3000, res = 400)
plotMD(qlf.UpdLgl_Time, main="LglIRDcr2+Upd 24H vs 96H", cex=0.75) + abline(h=c(-1, 1), col="blue")
dev.off()
tiff("Documents/NICDLgl_Time.tiff", height = 3000, width = 3000, res = 400)
plotMD(qlf.NICDLgl_Time, main="LglIRDcr2+NICD 24H vs 96H", cex=0.75) + abline(h=c(-1, 1), col="blue")
dev.off()






#DESEQ2 IS THE SHIT#

library(EnhancedVolcano)
library(airway)
library(magrittr)
library(gplots)
library(DESeq2)
library(lazyeval)
library(knitr)

targets<-read.csv("/home/dchatterjee/Documents/XFW_bulkRNASeqtumor_results/sample_info_tumor.csv")
group<-factor(paste(targets$Genotype,targets$Time,sep="."))
t <- cbind(targets,group=group)
x<-read.csv("/home/dchatterjee/Documents/XFW_bulkRNASeqtumor_results/tumor_gene_counts.csv", row.names="Geneid")
cat("Number of input rows:\n")
nrow(x)
x -> count.table
dds0<-DESeqDataSetFromMatrix(countData=x, colData=t, design = ~ group)
print(dds0)
slotNames(dds0)
dds.norm <-  estimateSizeFactors(dds0)
sizeFactors(dds.norm)

col.sample <- c("Lgl_RNAi.24H"="orange","Lgl_RNAi.96H"="red", "WT.24H"="gray77", "WT.96H"="gray29", "Lgl_RNAi_UAS_NICD.24H"="darkseagreen1", "Lgl_RNAi_UAS_NICD.96H"="darkseagreen4", "Lgl_RNAi_UAS_Upd.24H"="goldenrod", "Lgl_RNAi_UAS_Upd.96H"="goldenrod4", "UAS_Upd.24H"="firebrick1", "UAS_Upd.96H"="firebrick3", "UAS_NICD.24H"="royalblue1", "UAS_NICD.96H"="royalblue3")
t$color <- col.sample[as.vector(t$group)]
col.pheno.selected <- t$color

par(mfrow=c(3,1))

hist(as.matrix(count.table), col="blue", border="white", breaks=100)

hist(as.matrix(count.table), col="blue", border="white",
     breaks=20000, xlim=c(0,2000), main="Counts per gene",
     xlab="Counts (truncated axis)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(count.table + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

par(mfrow=c(1,1))

## Boxplots
expDesign <- t
boxplot(log2(count.table + epsilon), col=expDesign$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(Counts +1)")

## Density
## We will require one function from the affy package
if(!require("affy")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("affy")  
}
library(affy)
plotDensity(log2(count.table + epsilon), lty=1, col=expDesign$color, lwd=2)
grid()
legend("topright", legend=names(col.sample), col=col.sample, lwd=2)

prop.null <- apply(count.table, 2, function(x) 100*mean(x==0))
print(head(prop.null))
barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=expDesign$color, ylab='Samples', xlab='% of null counts')


par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)+epsilon),  col=col.pheno.selected, cex.axis=0.7,
las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)+epsilon),  col=col.pheno.selected, cex.axis=0.7,
las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts")
plotDensity(log2(counts(dds.norm)+epsilon),  col=col.pheno.selected,
xlab="log2(counts)", main="Density plot for Raw counts", cex.lab=0.7, panel.first=grid())
plotDensity(log2(counts(dds.norm, normalized=TRUE)+epsilon), col=col.pheno.selected,
xlab="log2(normalized counts)", main="Density plot for Normalized counts", cex.lab=0.7, panel.first=grid())

norm.counts <- counts(dds.norm, normalized=TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)
## sum(mean.counts==0) # Number of completely undetected genes
norm.counts.stats <- data.frame(
min=apply(norm.counts, 2, min),
mean=apply(norm.counts, 2, mean),
median=apply(norm.counts, 2, median),
max=apply(norm.counts, 2, max),
zeros=apply(norm.counts==0, 2, sum),
percent.zeros=100*apply(norm.counts==0, 2, sum)/nrow(norm.counts),
perc05=apply(norm.counts, 2, quantile, 0.05),
perc10=apply(norm.counts, 2, quantile, 0.10),
perc90=apply(norm.counts, 2, quantile, 0.90),
perc95=apply(norm.counts, 2, quantile, 0.95)
)
kable(norm.counts.stats)

par(mfrow=c(1,1))
mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
plot(x=log2(mean.counts), y=log2(variance.counts), pch=16, cex=0.5,
col=mean.var.col, main="Mean-variance relationship",
xlab="Mean log2(normalized counts) per gene",
ylab="Variance of log2(normalized counts)",
panel.first = grid())
abline(a=0, b=1, col="brown")

p <- 1/6 # the probability of success
n <- 1   # target for number of successful trials
# The density function
plot(0:10, dnbinom(0:10, n, p), type="h", col="blue", lwd=4)
dnbinom(0, n , p)
sum(dnbinom(0:5, n , p)) # == pnbinom(5, 1, p)
sum(dnbinom(0:10, n , p)) # == pnbinom(10, 1, p)
1-sum(dnbinom(0:10, n , p)) # == 1 - pnbinom(10, 1, p)
n <- 2
plot(0:30, dnbinom(0:30, n, p), type="h", col="blue", lwd=2,
main="Negative binomial density",
ylab="P(x; n,p)",
xlab=paste("x = number of failures before", n, "successes"))
# Expected value
q <- 1-p
(ev <- n*q/p)
abline(v=ev, col="darkgreen", lwd=2)
# Variance
(v <- n*q/p^2)
arrows(x0=ev-sqrt(v), y0 = 0.04, x1=ev+sqrt(v), y1=0.04, col="brown",lwd=2, code=3, , length=0.2, angle=20)
n <- 10
p <- 1/6
q <- 1-p
mu <- n*q/p
all(dnbinom(0:100, mu=mu, size=n) == dnbinom(0:100, size=n, prob=p))

## Performing estimation of dispersion parameter
dds.disp <- estimateDispersions(dds.norm)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
plotDispEsts(dds.disp)
alpha <- 0.0001

wald.test <- nbinomWaldTest(dds.disp)
res.DESeq2 <- results(wald.test, alpha=alpha, pAdjustMethod="BH")

## What is the object returned by nbinomTest()
class(res.DESeq2)
is(res.DESeq2) # a data.frame
slotNames(res.DESeq2)
head(res.DESeq2)
colnames(res.DESeq2)

par(mfrow=c(1,1))
hist(res.DESeq2$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(res.DESeq2$log2FoldChange, -log10(res.DESeq2$pvalue))
plot(res.DESeq2$log2FoldChange, -log10(res.DESeq2$padj), col=cols, panel.first=grid(),
main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res.DESeq2$log2FoldChange) > 2 & res.DESeq2$padj < alpha
text(res.DESeq2$log2FoldChange[gn.selected],
-log10(res.DESeq2$padj)[gn.selected],
lab=rownames(res.DESeq2)[gn.selected ], cex=0.4)

gn.most.sign <- rownames(res.DESeq2)[1]
gn.most.diff.val <- counts(dds.norm, normalized=T)[gn.most.sign,]
barplot(gn.most.diff.val, col=col.pheno.selected, main=gn.most.sign, las=2, cex.names=0.5)

gn.most.sign <- c("upd1")
gn.most.diff.val <- counts(dds.norm, normalized=T)[gn.most.sign,]
barplot(gn.most.diff.val, col=col.pheno.selected, main=gn.most.sign, las=2, cex.names=0.5)


## Perform the hierarchical clustering with
## A distance based on Pearson-correlation coefficient
## and average linkage clustering as agglomeration criteria
## We select gene names based on FDR (1%)
gene.kept <- rownames(res.DESeq2)[res.DESeq2$padj <= alpha & !is.na(res.DESeq2$padj)]
## We retrieve the normalized counts for gene of interest
count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
dim(count.table.kept)
if(!require("gplots")){
install.packages("gplots")
}
library("gplots")
## Perform the hierarchical clustering with
## A distance based on Pearson-correlation coefficient
## and average linkage clustering as agglomeration criteria
heatmap.2(as.matrix(count.table.kept),
scale="row",
hclust=function(x) hclust(x,method="average"),
distfun=function(x) as.dist((1-cor(t(x)))/2),
trace="none",
density="none",
labRow="",
cexCol=0.7)

dds0<-DESeqDataSetFromMatrix(countData=x, colData=t, design = ~ group)
dds <- DESeq(dds0, betaPrior=FALSE)

par(mfrow=c(2,2))
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

plotPCA(varianceStabilizingTransformation(dds), intgroup=c("Time")) + theme + theme(legend.position = "top")
plotPCA(varianceStabilizingTransformation(dds), intgroup=c("Genotype")) + theme + theme(legend.position = "top")
plotPCA(varianceStabilizingTransformation(dds), intgroup=c("gene_id")) + theme + theme(legend.position = "top")
plotPCA(varianceStabilizingTransformation(dds), intgroup=c("group")) + theme + theme(legend.position = "top")




#SELECT GENE COUNTS HEATMAP
select_var <- c('mre11', 'mei-9', 'mei-218', 'mei-41', 'mei-P26', 'mei-38', 'c(2)M', 'SMC2', 'SMC3', 'SMC5', 'mamo')
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
tiff("Documents/heatmap_select.tiff", height = 1200, width = 1700, res = 450)
pheatmap(highly_variable_lcpm, col=colorRampPalette(rev(brewer.pal(n = 10, name ="Spectral")))(100), kmeans_k=NA, fontsize=7, cluster_rows=TRUE, cluster_cols=TRUE,  clustering_distance_rows = "euclidean", cutree_cols=2)
dev.off()





#Enhanced Volcano Plots
#The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
dds0<-DESeqDataSetFromMatrix(countData=x, colData=t, design = ~ group)
dds <- DESeq(dds0, betaPrior=FALSE)

res1 <- results(dds, contrast = c('group', 'Lgl_RNAi.96H', 'Lgl_RNAi.24H'))
res1 <- lfcShrink(dds, contrast = c('group', 'Lgl_RNAi.96H', 'Lgl_RNAi.24H'), res=res1)
EnhancedVolcano(res1,
lab = rownames(res1),
title = 'Lgl-RNAi vs WT 24H',
x = 'log2FoldChange',
y = 'padj',
xlab = bquote(~Log[2]~ 'fold change'),
ylab = bquote(~-Log[10]~adjusted~italic(P)),
hline = c(10e-12, 10e-36, 10e-60, 10e-84),
hlineCol = c('grey0', 'grey25','grey50','grey75'),
hlineType = 'longdash',
hlineWidth = 0.8,
gridlines.major = FALSE,
gridlines.minor = FALSE,
xlim = c(-5, 8),
legend=c('NS','Log (base 2) fold-change','Adjusted p-value',
'Adjusted p-value & Log (base 2) fold-change'),
legendPosition = 'right',
legendLabSize = 16,
legendIconSize = 5.0)

EnhancedVolcano(res1,
                lab = rownames(res1),
                title = 'UAS-NICD vs w1118',
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('mus205', 'rad50', 'Top3alpha', 'mei-41', 'Xpc', 'spel1', 'lok', 'Msh6', 'Ku80', 'Gen', 'Blm', 'okr', 'CG32756', 'CG9272', 'Mlh1', 'mus308', 'mus101', 'RecQ5', 'Thd1', 'mei-9', 'Snm1', 'mre11', 'mus201', 'nbs', 'XRCC1', 'RpA-70', 'Fancm', 'Dcr-1', 'tefu', 'Uba2', 'Irbp', 'RfC4'),
                col=c('black', 'black', 'red1', 'red3'),
                colAlpha = 0.2,
                pCutoff = 0.05,
                FCcutoff = 1.8,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                hline = c(10e-12, 10e-36, 10e-60, 10e-84),
                hlineCol = c('grey0', 'grey25','grey50','grey75'),
                hlineType = 'longdash',
                hlineWidth = 0.8,
                boxedLabels = TRUE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                xlim = c(-5, 8),
                legend=c('NS','Log (base 2) fold-change','Adjusted p-value',
                         'Adjusted p-value & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30')
