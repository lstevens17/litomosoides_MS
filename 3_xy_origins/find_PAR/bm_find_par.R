library(ggplot2)
library(tidyverse)
library(patchwork)

# read in male cov and female cov
mcov <- read.table("bm_SRR3111510_cov_100k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
mtwon <- median(subset(mcov, chr != "Bm_v4_ChrX_scaffold_001")$cov)
mcov <- mcov %>% mutate(mtrans_cov = cov / mtwon)
mcov <- mcov %>% select(!(cov))
mcov <- mcov %>% filter(chr == "Bm_v4_ChrX_scaffold_001")

fcov <- read.table("bm_SRR488256_cov_100k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
ftwon <- median(subset(fcov, chr != "Bm_v4_ChrX_scaffold_001")$cov)
fcov <- fcov %>% mutate(ftrans_cov = cov / ftwon)
fcov <- fcov %>% select(!(cov))
fcov <- fcov %>% filter(chr == "Bm_v4_ChrX_scaffold_001")

# write bed file to find PAR
nonfemale_dominated_bed <- fcov %>% full_join(mcov) %>% mutate(ratio = mtrans_cov/ftrans_cov) %>% filter(ratio >= 0.8) %>% select((chr:stop))
write.table(nonfemale_dominated_bed, "bm_nonfemale_dominated.bed", sep="\t", quote = FALSE, col.names=F, row.names=F)
