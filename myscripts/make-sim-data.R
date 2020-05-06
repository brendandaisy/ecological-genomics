require(tidyverse)
require(DESeq2)
require(ggthemes)

n_gene <- 50
n_spt <- 10
cont_treat <- c(0, 1, 1.3, 2, 2.5, 3)

plot_contrast <- function(gene, discrete=FALSE, lfc=NULL) {
    if (discrete) {
        pred_points <- geom_point(aes(x=conds, y=effect), data=lfc, inherit.aes=FALSE, col='red', size=2, shape=2)
    }
    else {
        pred_points <- geom_abline(slope=res$log2FoldChange[gene], intercept=0, col='red')
    }
    tibble(conds, effect=log2_q[gene,]) %>%
        ggplot(aes(x=conds, y=effect)) +
        geom_point() +
        geom_line() +
        pred_points +
        theme_tufte()
}

K <- matrix(nrow=n_gene, ncol=n_spt * length(cont_treat))

log2_q <- matrix(rep(cont_treat, each=n_gene * n_spt), nrow=n_gene, ncol=n_spt * length(cont_treat))
disp <- rep(0, n_gene)
for (i in seq_len(n_gene)) {
    for (j in seq_len(n_spt * length(cont_treat))) {
        log2_q[i, j] <- if (i %% 2 == 0) log2_q[i, j]^2 else plogis(log2_q[i, j], 1.5, .5) * (i %% 5) * 2
    }
}

for (i in seq_len(n_gene)) {
    disp[i] <- rlnorm(1, meanlog=4 / mean(2^(log2_q[i,])) + 0.01, sdlog=0.5)
    for (j in seq_len(n_spt * length(cont_treat))) {
        K[i, j] <- rnbinom(1, mu=2^(log2_q[i, j]), size=disp[i])
    }
}

conds <- rep(cont_treat, each=n_spt)

### VISUALIZE THE DATA

tibble(conds, effect=log2_q[7,]) %>%
    ggplot(aes(x=conds, y=effect)) +
    geom_point() +
    geom_line() +
    theme_tufte()

tibble(g=rowMeans(K),
       d=disp,
       t=rep(c('thresh', 'sq'), times=length(disp) %/% 2)) %>%
    ggplot(aes(x=g, y=d, col=t)) +
    geom_point() +
    stat_function(fun=function(x) 4/x + 0.01, col='blue', linetype='dashed', n=300) +
    scale_x_continuous(trans='log2') +
    scale_y_continuous(trans='log2') +
    theme_tufte() +
    labs(x='Expected counts', y='Dispersion')

vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("conds"), returnData=FALSE) +
    theme_tufte()

qplot(conds, K[13,]) +
    theme_tufte() +
    labs(y='Counts for logistic response')

### RUN DESEQ2

dds <- DESeqDataSetFromMatrix(countData=K, colData=DataFrame(conds), design=~conds)
sizeFactors(dds) <- rep(1, 60)

dds <- DESeq(dds)

res <- results(dds)
resultsNames(dds)

plot_contrast(9)

### now seperate into discrete factor

conds <- as.factor(rep(as.character(cont_treat), each=n_spt))

dds <- DESeqDataSetFromMatrix(countData=K, colData=DataFrame(conds), design=~conds)
sizeFactors(dds) <- rep(1, 60)

dds <- DESeq(dds)
res <- results(dds)

r <- c(0)
for (s in resultsNames(dds)[2:6]) {
    l <- results(dds, name=s)$log2FoldChange[7]
    r <- c(r, l)
}

tibble(conds, effect=log2_q[7,]) %>%
    ggplot(aes(x=conds, y=effect)) +
    geom_point() +
    geom_line() +
    geom_point(aes(x=name, y=value), data=enframe(r), inherit.aes=FALSE, col='red', size=2, shape=2) +
    theme_tufte()
