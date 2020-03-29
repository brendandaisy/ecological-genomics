
--- title: "Brendan's Github Notebook"

output: html_document editor_options: chunk_output_type: console ---

## Author: Brendan Case
### Affiliation: UVM Computer Science
### E-mail contact: bcase@uvm.edu
### Website: bcase.w3.uvm.edu


### Start Date: 2020-01-13
### End Date: 2020-05-08
### Project Descriptions:


# Table of Contents:
* [Entry 1: 2020-01-13, Monday](#id-section1)
* [Entry 2: 2020-01-14, Tuesday](#id-section2) ------
<div id='id-section1'/>   

### Entry 1: 2020-01-22, Wednesday.



------ <div id='id-section2'/>

### Entry 2: 2020-01-29, Wednesday.

Useful markdown tips:

1. Embedding code: use backticks.

```
cd my data
ll
```

First 4 lines of a fastq file:

```
@GWNJ-0842:368:GW1809211440:2:1101:18355:2170 2:N:0:NCGTATCA+NTTCGCCT
NACAGAATTAATTTCATGGGAGAGCAAATGTACGTAACTAGGAAAATTTAAAGGGCTGTATGCACCCATATCCTCCTCATCCTCATAAGCATCATCAACAGATATGATAGGAAGAATTATGCCAGTTACAAGCTAGAACCCTCCCTTTTT
+
#AA<7A-A-AJA--<F--7F---7--<-<----<F----<-7--<<FFJ--<F-<--7--<F7-FJ<AF-A777<<F--<FF<-F--<-7-F-AFF<-F-<A-7-77-7A<F-7<-7<-77F-F<<F<-FFFJ-------)-)))77-AF
```

1st line contains metainfo on the sequence and the barcode on the
right

2nd line is basepairs (N = could be anything)

3rd line is a seperator

4th line is code for the "quality" of each basepair (i.e. was able to
read properly), corresponding to a p-value

Brendan: in charge of all the `PRK` files!!

------ <div id='id-section3'/>

### Entry 3: 2020-2-3 Monday.

Info updates with Kerry

Sewall Wright proposed the equation

$N_e = \frac{4 N_m N_f}{N_m + N_f}$

Another option is to estimate expected nucleotide diversity as $\pi =
4 N_e \mu$, where $\mu$ is the mutation rate per site per generation
(i.e. bit flip in an EA). (Recall that $\pi$ is defined over pairwise
SNPs $\pi_ij$ as $\pi = \sum_{ij} x_i x_j \pi_ij$, where $x_i$ is the
frequency of allele in $i$ in total, and $\pi_ij$ is no. SNPs between
2 individuals

Usually also interested in estimating $F_{ST}$. Maybe the simplest way
to do this is $F_{ST} = 1/(\pi +1)$. Another is $1 - H_s/H_t$, where
$H_s$ is within group heterozygosity, and $H_t$ is total het. in
population.

The site frequency spectrum is the distribution of the number of
individuals different from the dominant at each allelle (I.e. a
histogram of allelles, where each goes into the bin corresp. to no. of
individuals with the rare nucleotide).

Tajima's D attempts to use $pi$ to understand which alleles are under
natural selection.

<div id='id-section4'/>

### Entry 4: Feb 19th

Folding the SFS: for when the SFS is bimodal and you don't have enough
ancestry info to resolve whether your reference is ancestral.

There are other applications of SFS foldings as well

#### Estimating the (rough) SFS

Takes SNPs not as ground truth, but RVs with variance.

Starts with SFS as a prior, then calculate $\theta = 4 N_e \mu$ with
`doTheta` in bash, then use the `ThetaStats` function.

#### TODO: for next time

Update the git repo with all your scripts and results!

<div id='id-section5'/>

### Entry 5: Feb 26th

Notes from pipeline: we will be using `Salmon` to get alignment and
expression counts all at once (easier than getting sam files first),
and `DESeq2` for visualization and analysis (popular option; also, there
`EdgeR`.

When processing the Red spruce samples, seedlings were fully ground up
(roots through needles), so we can't pick up on expression at tissue
level.

We may remove day 5 from our statistics and only compare days 0 and
10, since 5 only got 1-4 reps per group.

3' tag sequencing requires very few samples (though need high quality
RNA). Our lab here is expensive and has some other issues, so common
to use Cornell's service here.

Part of what makes RNA more volatile is the O-H end on RNA not found
on DNA.

Fractors/num. levels:

- Treatment (3)
- SourceClim (2)
- Time (3)

Questions:

- Do ind from diff climsources have diff gene exp
- Expression model; exp ~ time + SC + treatment + time x SC + ... + (family as random variable?)

### Entry 6

Starting with a data table for doing DESeq. The table consists of the
number of reads which mapped to a certain gene (row), by which fastq
file the read came from (col).


