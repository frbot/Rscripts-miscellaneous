## DESeq2 Analysis based on DESeq2 vignette

## Load libraries
library("airway")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("vsn")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("genefilter")
library("tibble")


## Set wd containing sample_table.csv and bam files
indir = (".")
csvfile <- file.path(indir, "sample_table.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
filenames <- file.path(indir, paste0(sampleTable$Sample, "-sorted.bam"))

## Handle BAM files
bamfiles <- BamFileList(filenames, yieldSize=2000000)

## Defining gene model
gtffile <- file.path(indir,"Genemodel.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
ebg <- cdsBy(txdb, by="gene")

## Read counting step
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=TRUE)

colData(se) <- DataFrame(sampleTable)

## Run DESeq2
se$Condition <- relevel(se$Condition, "Control")
dds <- DESeqDataSet(se, design = ~ Condition)

## Prefiltering the dataset removing the columns with 0
dds <- dds[ rowSums(counts(dds)) > 1, ]

## Poisson counts  with a range of lambda from 0.1 to 100
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
png("Counts-before-transform.png", 1000, 1000, pointsize=20)
meanSdPlot(cts, ranks = FALSE)
dev.off()

## Distribution for log transformed counts
log.cts.one <- log2(cts + 1)

#If you want to use VST transforation
#vsd <- vst(dds, blind = FALSE)
#rld <- rlog(dds, blind = FALSE)

## Show the effect of transformation on the data
rld <- rlog(dds, blind = FALSE)
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")
png("Counts-log-transformed.png", 1000, 1000, pointsize=20)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
dev.off()

## Compute sample to sample distances and plot the result
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Condition, rld$Experiment, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

## PCA plot
#plotPCA(rld, intgroup = c("Condition", "Experiment"))

pcaData <- plotPCA(rld, intgroup = c( "Condition", "Experiment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("Plot-PCA.png", 1000, 1000, pointsize=20)
ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition, shape = Experiment)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()

## Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)

## Load annotation
annotation <- file.path(indir, "Annotation.csv")
annot <- read.csv(annotation)

## DFR < 0.05 filtered results
res2 <- subset(res, res$padj < 0.05)
res2 <- res2[order(res2$padj), ]
table(res2$padj<0.05)

resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
resannot2 <- resdata2 %>% left_join(annot)
write.csv(resannot2, file="diffexpr-results-filtered.csv")

## Write full results with linked annotation
sum(res$padj<0.05)
table(res$padj<0.05)
res <- res[order(res$padj), ]

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
resannot <- resdata %>% left_join(annot)
write.csv(resannot, file="diffexpr-results-full.csv")

## Examine plot of p-values
png("Pvalues-hist.png", 1000, 1000, pointsize=20)
hist(res$pvalue, breaks=50, col="grey")
dev.off()

#Dispersion estimate
png("Dispersion.png")
plotDispEsts(dds, ylim = c(1e-6, 1e1) )
dev.off()

#MA-plot
res <- lfcShrink(dds, coef="Condition_Test_vs_Control", type="apeglm")
png("MAplot.png", 1000, 1000, pointsize=20)
plotMA(res, ylim = c(-5, 5))
dev.off()

#Gene Clustering

mat = assay(rld)[ head(order(res$padj),50), ]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld)[,c("Condition")])
colnames(df) = "Condition"
rownames(df) = colnames(mat)
pdf("Heatmap-top50.pdf")
pheatmap(mat, annotation_col=df, fontsize_row = 7)
dev.off()

png("Heatmap-top50.png", 1000, 1000, pointsize=20)
pheatmap(mat, annotation_col=df, fontsize_row = 7)
dev.off()

mat = assay(rld)[ head(order(res$padj),100), ]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld)[,c("Condition")])
colnames(df) = "Condition"
rownames(df) = colnames(mat)
pdf("Heatmap-top100.pdf")
pheatmap(mat, annotation_col=df, fontsize_row = 5)
dev.off()

png("Heatmap-top100.png", 1000, 1000, pointsize=20)
pheatmap(mat, annotation_col=df, fontsize_row = 5)
dev.off()

mat = assay(rld)[ head(order(res$padj),200), ]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld)[,c("Condition")])
colnames(df) = "Condition"
rownames(df) = colnames(mat)
pdf("Heatmap-top200.pdf")
pheatmap(mat, annotation_col=df, fontsize_row = 5)
dev.off()

png("Heatmap-top200.png", 1000, 1000, pointsize=20)
pheatmap(mat, annotation_col=df, fontsize_row = 5)
dev.off()
