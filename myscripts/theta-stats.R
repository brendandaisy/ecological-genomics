require(tidyverse)

setwd('c:/Users/brendandaisy/Documents/phd/ecological-genomics/myresults/')

list.files()

sfs <- scan('AB_outFold.sfs') # PRK is my pop for the class

(pct_poly <- 100*(1 - (sfs[1] / sum(sfs))))

div <- read.table('AB_folded_allsites.thetas.idx.pestPG') %>%
    as_tibble

colnames(div) <- c('window', 'chrname', 'wincenter', 'tW', 'tP', 'tF', 'tH', 'tL', 'tajD', 'fulif', 'fuliD', 'fayH', 'zengsE', 'numSites')

div$tWpersite <- div$tW/div$numSites
div$tPpersite <- div$tP/div$numSites

hist(div$tWpersite)
hist(div$tPpersite)
hist(div$tajD)
barplot(sfs[-1])

summary(div)
     
