#!/bin/bash

# set the repo

myrepo="/users/b/c/bcase/ecological-genomics"
mypop="PRK"
outfolder="${myrepo}/myresults/bamstats/"

for file in ${myrepo}/myresults/BWA/${mypop}*sorted.rmdup.bam
do
    f=${file/.sorted.rmdup.bam/} # removes the stuff between /./
    name=`basename ${f}`
    echo ${name} >> ${outfolder}/${mypop}.names.txt
    # 1. calculate stats on the mapping flags
    echo "Num.reads R1 R2 Paired MateMapped Singletons MateMappedDiffChr" > ${outfolder}${mypop}.flagstats.txt
    samtools flagstat ${file} | awk 'NR>=6&&NR<=12 {print $1}' | column -x >> ${outfolder}/${mypop}.flagstats.txt

    # 2. calculate depth of coverage from bam files (important for statistics on reads)
    samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}' > ${outfolder}/${mypop}.coverage.txt
done
