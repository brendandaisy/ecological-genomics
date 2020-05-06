require(DESeq2)
require(tidyverse)

## Import the counts matrix made using Salmon
counts <- read.table("../mydata/RS_counts_samples/RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1) %>%
    round # Need to round because DESeq wants only integers

## Import the samples description table - links each sample to factors of the experimental design.
## Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("../mydata/RS_counts_samples/RS_samples.txt",
                    header=TRUE,
                    stringsAsFactors = TRUE,
                    row.names=1,
                    colClasses=c('factor', 'factor', 'factor', 'factor'))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = conds, 
                              design = ~day + treatment)

## get the sample-wise size factors (each gene within a sample has same size factors)
dd_sf <- estimateSizeFactors(dds)

## get gene-wise dispersion parameters
dd_disp <- estimateDispersions(dd_sf)
dispersions(dd_disp) # the NAs corresp to counts per gene are 0 (seems to depend on the model matrix!)
