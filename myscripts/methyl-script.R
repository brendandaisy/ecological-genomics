require(methylKit)
require(tidyverse)
require(pheatmap)

# first, we want to read in the raw methylation calls with methylkit

## set directory with absolute path (why is this necessary? I have no idea, but gz files wont work with relative paths)

dir <- 'c:/Users/brendandaisy/Documents/phd/ecological-genomics/mydata/Epigenetics_data'

## read in the sample ids
samples <- read.table(paste0(dir, '/sample_id.txt'),
                      header=FALSE)

## now point to coverage files
files <- file.path(dir, samples$V1)

## convert to list

## get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz","",samples$V1))

## use methRead to read in the coverage files
myobj <- methRead(location=as.list(files),
                  sample.id = nmlist,
                  assembly='atonsa',
                  dbtype='tabix',
                  context='CpG',
                  resolution='base',
                  mincov=20,
                  treatment=rep(0:4, each=4),
                  pipeline='bismarkCoverage',
                  dbdir='c:/Users/brendandaisy/Documents/phd/ecological-genomics/myresults')

######
## visualize coverage and filter
######

## We can look at the coverage for individual samples with getCoverageStats()
getCoverageStats(myobj[[1]], plot=TRUE)

## and can plot all of our samples at once to compare.

## filter samples by depth with filterByCoverage()
filtered.myobj <- filterByCoverage(myobj,
                                   lo.count=20,
                                   lo.perc=NULL,
                                   hi.count=NULL,
                                   hi.perc=97.5)

######
## merge samples
######

##Note! This takes a while and we're skipping it

## use unite() to merge all the samples. We will require sites to be present in each sample or else will drop it

meth <- methylKit::unite(filtered.myobj,
                         mc.cores=3,
                         suffix='united',
                         db.dir='c:/Users/brendandaisy/Documents/phd/ecological-genomics/myresults')

meth <- methylKit:::readMethylBaseDB(
     dbpath = paste0(dir, "/methylBase_united.txt.bgz"),
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

                                        # percMethylation() calculates the percent methylation for each site and sample

pm <- percMethylation(meth)

#plot methylation histograms
ggplot(gather(as.data.frame(pm)), aes(value)) + 
    geom_histogram(bins = 10, color="black", fill="grey") + 
    facet_wrap(~key, ncol=4)

                                        # calculate and plot mean methylation
sp.means <- colMeans(pm)
p.df <- data.frame(sample=names(sp.means),
          group = substr(names(sp.means), 1,6),
          methylation = sp.means)

ggplot(p.df, aes(x=group, y=methylation, color=group)) + 
    stat_summary(color="black") +
    geom_jitter(width=0.1, size=3) 

# sample clustering
clusterSamples(meth,
               dist='correlation',
               method='ward.D',
               plot=TRUE)

# PCA
PCASamples()

                                        # subset with reorganize()


meth_sub <- reorganize(meth,
                sample.ids =c("AA_F00_1","AA_F00_2","AA_F00_3", "AA_F00_4",
                              "HH_F25_1","HH_F25_2","HH_F25_3","HH_F25_4"),
                treatment = c(0,0,0,0,1,1,1,1),
                save.db=FALSE)
                             
# calculate differential methylation

myDiff <- calculateDiffMeth(meth_sub,
                            overdispersion='MN',
                            mc.cores=1,
                            suffix='AA_HH',
                            adjust='qvalue',
                            test='Chisq')

                                        # get all differentially methylated bases

myDiff <- getMethylDiff(myDiff,
                        qvalue=.05,
                        difference=10)

                                        # we can visualize the changes in methylation frequencies quickly.

hist(getData(myDiff)$meth.diff)

# get hyper methylated bases

                                        # get hypo methylated bases

#heatmaps first

# get percent methylation matrix

# make a dataframe with snp id's, methylation, etc.

# add snp, chr, start, stop


####
# heatmap
####

pm <- percMethylation(meth_sub)
# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in,]

# add snp, chr, start, stop

din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, getData(myDiff)$start, sep=":"), din, pm.sig)

colnames(df.out) <- c("snp", colnames(din), colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[,5:7]))

pheatmap(pm.sig, show_rownames=FALSE)

                                        # we can also normalize

ctrmean <- rowMeans(pm.sig[,1:4])
h.norm <- pm.sig - ctrmean
pheatmap(h.norm, show_rownames=FALSE)

#####
#let's look at methylation of specific snps
####

# convert data frame to long form

## write bed file for intersection with genome annotation

