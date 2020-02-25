#!/bin/bash

# set the repo and output folder

myrepo="/users/b/c/bcase/ecological-genomics"
mkdir ${myrepo}/myresults/angsd
outfolder="${myrepo}/myresults/angsd"
mypop="PRK"

# put the bam file names in the angsd folder 
ls /data/project_data/RS_ExomeSeq/mapping/BWA/${mypop}*sorted.rm*.bam > ${outfolder}/${mypop}_bam.list

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# Estimating GL's and allele frequencies for all sites with ANGSD

ANGSD -b ${outfolder}/${mypop}_bam.list \
      -ref ${REF} -anc ${REF} \
      -out ${outfolder}/${mypop}_folded_allsites \
      -nThreads 1 \
      -remove_bads 1 \
      -C 50 \
      -baq 1 \
      -minMapQ 20 \
      -minQ 20 \
      -setMinDepth 3 \
      -minInd 2 \
      -setMinDepthInd 1 \
      -setMaxDepthInd 17 \
      -skipTriallelic 1 \
      -GL 1 \
      -doCounts 1 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -doSaf 1 \
      -doHWE 1
# -SNP_pval 1e-6

# now do it with folding since we have been told a priori the SFS was bimodal

ANGSD -b ${outfolder}/${mypop}_bam.list \
      -ref ${REF} -anc ${REF} \
      -out ${outfolder}/${mypop}_folded_allsites \
      -nThreads 1 \
      -remove_bads 1 \
      -C 50 \
      -baq 1 \
      -minMapQ 20 \
      -minQ 20 \
      -setMinDepth 3 \
      -minInd 2 \
      -setMinDepthInd 1 \
      -setMaxDepthInd 17 \
      -skipTriallelic 1 \
      -GL 1 \
      -doCounts 1 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -doSaf 1 \
      -fold 1 # this will fold the SFS

# Nwo to get a rough first estimate of the SFS, use org SFS as prior

realSFS ${outfolder}/${mypop}_folded_allsites.saf.idx -maxIter 1000 -tole 1e-6 -P 1 > ${outfolder}/${mypop}_outFold.sfs

# get refined est of SFS and doTheta

ANGSD -b ${outfolder}/${mypop}_bam.list \
      -ref ${REF} -anc ${REF} \
      -out ${outfolder}/${mypop}_folded_allsites \
      -nThreads 1 \
      -remove_bads 1 \
      -C 50 \
      -baq 1 \
      -minMapQ 20 \
      -minQ 20 \
      -setMinDepth 3 \
      -minInd 2 \
      -setMinDepthInd 1 \
      -setMaxDepthInd 17 \
      -skipTriallelic 1 \
      -GL 1 \
      -doCounts 1 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -doSaf 1 \
      -fold 1 \
      -pest ${outfolder}/${mypop}_outFold.sfs \
      -doThetas 1

# use the doTheta out from above to est nucleotide diversity

thetaStat do_stat ${outfolder}/${mypop}_folded_allsites.thetas.idx
      
