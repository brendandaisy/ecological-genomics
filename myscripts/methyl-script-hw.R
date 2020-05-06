require(methylKit)
require(tidyverse)
require(pheatmap)
require(ggthemes)

## set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)

dir <- 'c:/Users/brendandaisy/Documents/phd/ecological-genomics/mydata/Epigenetics_data'

## load united samples (remember to specify how filtering was done before uniting)
meth <- methylKit:::readMethylBaseDB(dbpath = paste0(dir, "/methylBase_united.txt.bgz"),
                                     dbtype = "tabix",
                                     sample.id =   unlist(nmlist),
                                     assembly = "atonsa", # this is just a string. no actual database
                                     context = "CpG",
                                     resolution = "base",
                                     treatment = c(0,0,0,0,
                                                   1,1,1,1,
                                                   2,2,2,2,
                                                   3,3,3,3,
                                                   4,4,4,4),
                                     destrand = FALSE)

## get avg methyation in each treatment group
percMethylation(meth) %>%
    as_tibble %>%
    gather(label, pm) %>%
    mutate(treat = str_sub(label, 1, -3)) %>%
    group_by(treat) %>%
    summarize(mean_pm = mean(pm))

clusterSamples(meth)

meth_gen_cont <- reorganize(meth,
                            sample.ids =c("AA_F00_1","AA_F00_2","AA_F00_3", "AA_F00_4",
                                          "HH_F25_1","HH_F25_2","HH_F25_3","HH_F25_4"),
                            treatment = c(0,0,0,0,1,1,1,1),
                            save.db=FALSE)

PCASamples(meth_gen_cont)

meth_cont <- reorganize(meth,
                        sample.ids = c("AA_F25_1","AA_F25_2", "AA_F25_3", "AA_F25_4",
                                       "AH_F25_1","AH_F25_2", "AH_F25_3", "AH_F25_4"),
                        treatment = rep(c(0, 1), each=4),
                        save.db = FALSE)

PCASamples(meth_cont)

diff_meth <- calculateDiffMeth(meth_cont,
                               overdispersion='MN',
                               mc.cores=1,
                               suffix='AA_AH',
                               adjust='qvalue',
                               test='Chisq')

## get all differentially methylated bases

diff_base <- getMethylDiff(diff_meth,
                           qvalue=.05,
                           difference=10)

hist(getData(diff_base)$meth.diff)

getMethylDiff(diff_meth,
              qvalue=.05,
              difference=10,
              type='hyper') %>%
    nrow

getMethylDiff(diff_meth,
              qvalue=.05,
              difference=10,
              type='hypo') %>%
    nrow

## make a dataframe with snp id's, methylation, etc.
pm <- percMethylation(meth_cont)
sig.in <- as.numeric(row.names(diff_base))
pm.sig <- pm[sig.in,]

## add snp, chr, start, stop
din <- getData(diff_base)[,1:3]
dfb_out <- cbind(paste(getData(diff_base)$chr, getData(diff_base)$start, sep=":"), din, pm.sig)
colnames(dfb_out) <- c("snp", colnames(din), colnames(dfb_out[5:ncol(dfb_out)]))
dfb_out <- (cbind(dfb_out, getData(diff_base)[,5:7]))

## heatmap
ctrmean <- rowMeans(pm.sig[,1:4])
h.norm <- (pm.sig-ctrmean)
pheatmap(h.norm, show_rownames = FALSE)

## look at ind. SNPs
dfb_out
df.plot <- dfb_out[,c(1,5:12)] %>% pivot_longer(-snp, values_to = "methylation")
df.plot$group <- substr(df.plot$name,1,2)
head(df.plot)

means <- filter(df.plot, snp %in% c("LS274951.1:141", 'LS328870.1:810', 'LS068114.1:6322')) %>%
    mutate(group = paste(group, '_F25')) %>%
    group_by(group, snp) %>%
    summarize(mean = mean(methylation)) %>%
    ungroup

filter(df.plot, snp %in% c("LS274951.1:141", 'LS328870.1:810', 'LS068114.1:6322')) %>%
    mutate(group = paste(group, '_F25')) %>%
    ggplot(aes(x=group, y=methylation, color=snp)) +
    ## stat_summary(fun.data = "mean_se", size = 1.8) +
    geom_jitter(size=2, width=.1) +
    geom_line(aes(y = mean, group=snp), data=means) +
    geom_point(aes(y = mean, group=snp, fill=snp), data=means, size=3.2, shape=21, col='black') +
    labs(x='Treatment', y='% Methylation', title='Methylation rate for enriched SNPs', col='SNP location', fill='SNP location') +
    theme_tufte()

ggplot(means, aes(x=group, y=mean, group=snp, col=snp)) +
    geom_point() +
    geom_line()

## write bed file for intersection with genome annotation
write.table(file = "../myresults/diffmeth-aa-ah.bed",
            data.frame(chr= dfb_out$chr, start = dfb_out$start, end = dfb_out$end),
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep="\t")

read.table('../myresults/hits.bed')
