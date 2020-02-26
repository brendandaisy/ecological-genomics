setwd('c:/Users/brendandaisy/Documents/phd/ecological-genomics/myresults/angsd')

list.files()

sfs <- scan('PRK_outFold.sfs')

(total_sites <- sum(sfs)) # get the total number of sites
(pct_poly <- 100*(1 - (sfs[1] / total_sites))) # calculate percentage of sites which were SNPs

div <- read.table('PRK_folded_allsites.thetas.idx.pestPG')

colnames(div) <- c('window', 'chrname', 'wincenter', 'tW', 'tP', 'tF', 'tH', 'tL', 'tajD', 'fulif', 'fuliD', 'fayH', 'zengsE', 'numSites')

div$tWpersite <- div$tW/div$numSites
div$tPpersite <- div$tP/div$numSites

par(mfrow=c(2, 2))
hist(div$tWpersite, main='Estimated Theta (per site)')
hist(div$tPpersite, main='Observed Pi (per site)')
hist(div$tajD, main='Tajima\'s D (per window?)')
barplot(sfs[-1], main='Site Frequency Spectrum')

summary(div)
