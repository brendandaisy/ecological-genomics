library(DESeq2)
require(tidyverse)
library(scales)
library(ggpubr)
library(wesanderson)
require(ggthemes)
​
## Import the counts matrix made using Salmon
countsTable <- read.table("../mydata/RS_counts_samples/RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers
head(countsTableRound)
​
## Import the samples description table - links each sample to factors of the experimental design.
## Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("../mydata/RS_counts_samples/RS_samples.txt",
                    header=TRUE,
                    stringsAsFactors = TRUE,
                    row.names=1,
                    colClasses=c('factor', 'factor', 'factor', 'factor'))
​
### Try with only Day 10 data
​
day10countstable <- countsTableRound %>%
    select(contains("10"))
​
conds10 <- subset(conds, day=="10")
​
## Let's see how many reads we have from each sample:
barplot(colSums(countsTableRound), las=3, cex.names=0.5,names.arg = substring(colnames(countsTableRound),1,13))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd =2)
​
## What's the average number of counts per gene
mean(rowSums(day10countstable))
median(rowSums(day10countstable))
## wow! This shows dispersion across genes - differences in magnitude of expression

​
## Create a DESeq object and define the experimental design
dds <- DESeqDataSetFromMatrix(countData = day10countstable,
                              colData = conds10, 
                              design = ~ climate + treatment + climate:treatment)

# Filter out genes with few reads 
dds <- dds[rowSums(counts(dds)) > 76]

dds <- DESeq(dds)​

## List the results you've generated
resultsNames(dds)

## find avgs after normalization
mean(rowSums(counts(dds, normalized=TRUE)))
median(rowSums(counts(dds, normalized=TRUE)))

## Order and list and summarize results from specific contrasts. Here
## you set your adjusted p-value cutoff, can make summary tables of
## the number of genes differentially expressed (up- or
## down-regulated) for each contrast

res_climate <- results(dds, alpha = 0.05, name='climate_HD_vs_CW') %>%
    as_tibble(rownames='gene') %>%
    filter(pvalue <= 0.05) %>%
    arrange(pvalue)

res_treat_d <- results(dds, alpha = 0.05, name='treatment_D_vs_C') %>%
    as_tibble(rownames='gene') %>%
    filter(pvalue <= 0.05) %>%
    arrange(pvalue)

res_treat_h <- results(dds, alpha = 0.05, name='treatment_H_vs_C') %>%
    as_tibble(rownames='gene') %>%
    filter(pvalue <= 0.05) %>%
    arrange(pvalue)

## find all genes which had signif G or E effect
res_all_primary <- reduce(list(res_climate, res_treat_d, res_treat_h),
                              ~distinct(union(.x, .y), gene, .keep_all=TRUE))

res_d_int <- results(dds, alpha = 0.05, name='climateHD.treatmentD') %>%
    as_tibble(rownames='gene') %>%
    filter(pvalue <= 0.05) %>%
    arrange(pvalue)

res_h_int <- results(dds, alpha = 0.05, name='climateHD.treatmentH') %>%
    as_tibble(rownames='gene') %>%
    filter(pvalue <= 0.05) %>%
    arrange(pvalue)

res_all_int <- union(res_d_int, res_h_int) %>%
    distinct(gene, .keep_all=TRUE)
​
r1 <- setdiff(res_all_int$gene, res_all_primary$gene)

length(setdiff(r1, res_d_int$gene)) # those only in the h interaction
length(setdiff(r1, res_h_int$gene))
​
# PCA
vsd <- vst(dds, blind=FALSE)
​
pca_dat <- plotPCA(vsd, intgroup=c("climate","treatment"), returnData=TRUE)
vpercentVar <- round(100 * attr(pca_dat, "percentVar"))
​
pca_dat$treatment <- factor(pca_dat$treatment,
                         levels = c("C","H","D"),
                         labels = c("C","H","D"))
​
ggplot(pca_dat, aes(PC1, PC2, color=treatment, shape=climate)) +
    geom_point(size=3) +
    labs(x=paste0('PC1: ', vpercentVar[1], '%'),
         y=paste0('PC2: ', vpercentVar[2], '%'),
         color='Treatment',
         shape='Source climate') +
    theme_tufte()
​
## # Counts of specific top gene! (important validatition that the normalization, model is working)
## (d <- plotCounts(dds, gene="MA_7017g0010", intgroup = (c("treatment","climate","day")), returnData=TRUE))
## ​
## p <- ggplot(d, aes(x=treatment, y=count, color=day, shape=climate)) + 
##   theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
## p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
##   scale_x_discrete(limits=c("C","H","D"))
## p
## ​
## p <- ggplot(d, aes(x=treatment, y=count, shape=climate)) + 
##     theme_minimal() +
##     theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
## p
​
# Heatmap of top 20 genes sorted by pvalue
​
require(pheatmap)
topgenes <- head(rownames(res),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate")])
pheatmap(mat, annotation_col=df)
