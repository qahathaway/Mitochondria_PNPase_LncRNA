####SummarizedExperiment####
##Workflow for human and murine long non-coding RNA (lncRNA) sequencing in the
##mitochondria and cytoplasm

library("path_to_your_package")

indir <- system.file("extdata", package="path_to_your_package", mustWork=TRUE)

list.files(indir)

csvfile <- file.path(indir, "path_to_your_sample_table")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable

filenames <- file.path(indir, paste0(sampleTable$Run, ".bam"))
file.exists(filenames)

library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])

library("GenomicFeatures")

gtffile <- file.path(indir,"path_to_your_GTF_file")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

##Add Chromosome Initials

#seqlevels(txdb) <- sub("", "chr", seqlevels(txdb))
#seqlevels(txdb) <- sub("chrMT", "chrM", seqlevels(txdb))
#seqlevels(txdb)

ebg <- exonsBy(txdb, by="gene")
ebg

library("GenomicAlignments")
library("BiocParallel")

register(SerialParam())

##Either singleEnd=TRUE & fragments=FALSE (SE) or singleEND=FALSE & fragments=TRUE (PE)
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

se
dim(se)
assayNames(se)
head(assay(se), 5)
colSums(assay(se))
rowRanges(se)
str(metadata(rowRanges(se)))

####Make Counts Table####

#assay(se)
resOrderedCounts <- as.data.frame(assay(se))
write.csv(resOrderedCounts, file = "path_to_your_counts.csv")

colData(se)
colData(se) <- DataFrame(sampleTable)
colData(se)

library("magrittr")
se$dex %<>% relevel("untrt2")
se$dex
se$dex <- relevel(se$dex, "untrt2")


####Starting from SummarizedExperiment####

round( colSums(assay(se)) / 1e6, 1 )
colData(se)

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ dex)


####Starting from Count Matrices####

#HumanLncCounts_noncoding <- read.csv("path_to_your_counts.csv", row.names = 1)

#countdata <- HumanLncCounts_noncoding
#head(countdata, 3)

#coldata <- DataFrame(read.csv(file = "path_to_your_sample_table.csv"
#                              , row.names = 1))

#dds <- DESeqDataSetFromMatrix(countData = countdata,
#                                 colData = coldata,
#                                 design = ~ dex)



####Sample Distribution####

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 2, ]
nrow(dds)

##Data Normalization Assessment
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

##Sample Distances
sampleDists <- dist(t(assay(rld)))
sampleDists

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Run, rld$SampleName, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$Run, rld$SampleName, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

##PCA Plot
plotPCA(rld, intgroup = c("dex", "Run"))

pcaData <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

##MDS Plot
mds <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()


####Differential Expression Analysis####
dds <- DESeq(dds)

res <- results(dds)
res
summary(res)

##Middle group is "treatment", end group is "control"
res <- results(dds, contrast=c("dex","trt2","untrt2"), alpha = 0.05)
mcols(res, use.names = TRUE)
summary(res)

####Figures####


##VennDiagram
##Four groups, unique genes

library(limma)

res1 <- results(dds, contrast=c("dex","trt","untrt"), alpha = 0.05)
res1 <- res1[which(res1$padj < 0.05),]
res1.genes <- row.names(res1)

res2 <- results(dds, contrast=c("dex","trt2","untrt2"), alpha = 0.05)
res2 <- res2[which(res2$padj < 0.05),]
res2.genes <- row.names(res2)

res3 <- results(dds, contrast=c("dex","trt2","trt"), alpha = 0.05)
res3 <- res3[which(res3$padj < 0.05),]
res3.genes <- row.names(res3)

res4 <- results(dds, contrast=c("dex","untrt2","untrt"), alpha = 0.05)
res4 <- res4[which(res4$padj < 0.05),]
res4.genes <- row.names(res4)

Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes, res4.genes)))

res1.genes.2 <- Unique %in% res1.genes
res2.genes.2 <- Unique %in% res2.genes
res3.genes.2 <- Unique %in% res3.genes
res4.genes.2 <- Unique %in% res4.genes

counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2,res4.genes.2)

results.1 <- vennCounts(counts.1, include="both")

vennDiagram(results.1, include="both", cex = 1, names = c("","","",""), circle.col = c("blue", "red", "green", "black"))


##Other Visualizations

library(vidger)

vsBoxPlot(dds, d.factor = "dex", type = "deseq")
vsDEGMatrix(dds, d.factor = "dex", padj = 0.05, type = "deseq")
vsFourWay("trt2", "trt", "untrt2", d.factor = "dex", dds, type = "deseq")
vsMAMatrix(dds, d.factor = "dex", type = "deseq", y.lim = c(-10,10))
vsScatterMatrix(dds, d.factor = "dex", type = "deseq")
vsScatterPlot("untrt2", "trt2", dds, d.factor = "dex", type = "deseq")
vsVolcano("untrt2", "trt2", dds, d.factor = "dex", type = "deseq", x.lim = c(-10,10), padj = 0.05)
vsVolcanoMatrix(dds, d.factor = "dex", type = "deseq", lfc = 2, padj = 0.05, x.lim = c(-8,8),
                title = FALSE, legend = TRUE, grid = TRUE, counts = FALSE, facet.title.size = 10)

##Volcano Plots

library(EnhancedVolcano)

EnhancedVolcano(res2, lab = rownames(res), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-5, 7),  ylim = c(0, -log10(10e-11)))

##Strict Parameters
#res.05 <- results(dds, alpha = 0.01)
#table(res.01$padj < 0.01)
#resLFC1 <- results(dds, lfcThreshold=1)
#table(resLFC1$padj < 0.1)

#Plotting Results
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))

GeneofChoice <- "Ensembl_ID_gene_of_choice"
plotCounts(dds, gene = GeneofChoice, intgroup=c("dex"))

library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = GeneofChoice, intgroup = c("dex","Run"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = Run)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

##MA-Plot
res <- lfcShrink(dds, contrast=c("dex","trt2","untrt2"), res=res)
plotMA(res, ylim = c(-6, 6))

##MA-Plot with top gene
plotMA(res, ylim = c(-6, 6))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

##Histogram of p vaules for genes with mean normalized count larger than 2
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

##Gene Clustering

library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 500)
#or#
topGenes <- head(order(res$padj),decreasing = TRUE, 500)

mat  <- assay(rld)[ topGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)

####Annotating and Exporting####

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$padj),]
resOrdered <- res[order(res$baseMean, decreasing = TRUE),]
head(resOrdered, 8)

resOrderedDF <- as.data.frame(resOrdered)[1:60000, ]
write.csv(resOrderedDF, file = "path_to_file_differential_expression.csv")
