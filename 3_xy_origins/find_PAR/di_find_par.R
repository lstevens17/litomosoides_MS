library(ggplot2)
library(tidyverse)
library(patchwork)

# read in male cov and female cov
mcov <- read.table("di_SRR13154015_cov_100k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
mtwon <- median(subset(mcov, chr != "dirofilaria_immitis_chrX")$cov)
mcov <- mcov %>% mutate(mtrans_cov = cov / mtwon)
mcov <- mcov %>% select(!(cov))
mcov <- mcov %>% filter(chr == "dirofilaria_immitis_chrX")

fcov <- read.table("di_ERR034941_cov_100k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
ftwon <- median(subset(fcov, chr != "dirofilaria_immitis_chrX")$cov)
fcov <- fcov %>% mutate(ftrans_cov = cov / ftwon)
fcov <- fcov %>% select(!(cov))
fcov <- fcov %>% filter(chr == "dirofilaria_immitis_chrX")

# write bed file to find PAR
nonfemale_dominated_bed <- fcov %>% full_join(mcov) %>% mutate(ratio = mtrans_cov/ftrans_cov) %>% filter(ratio > 0.8) %>% select((chr:stop))
write.table(nonfemale_dominated_bed, "di_nonfemale_dominated.bed", sep="\t", quote = FALSE, col.names=F, row.names=F)
