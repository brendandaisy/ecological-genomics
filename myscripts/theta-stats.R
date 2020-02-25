setwd('c:/Users/brendandaisy/Documents/phd/ecological-genomics/myresults/')

list.files()

sfs <- scan('PRK_outFold.sfs')

(total_sites <- sum(sfs)) # get the total number of sites
(pct_poly <- 100*(1 - (sfs[1] / total_sites))) # calculate percentage of sites which were SNPs

div <- read.table('PRK_folded_allsites.thetas.idx.pestPG')

colnames(div) <- c('window', 'chrname', 'wincenter', 'tW', 'tP', 'tF', 'tH', 'tL', 'tajD', 'fulif', 'fuliD', 'fayH', 'zengsE', 'numSites')

div$tWpersite <- div$tW/div$numSites
div$tPpersite <- div$tP/div$numSites

hist(div$tWpersite)
hist(div$tPpersite)
hist(div$tajD)
barplot(sfs[-1])

summary(div)
     
