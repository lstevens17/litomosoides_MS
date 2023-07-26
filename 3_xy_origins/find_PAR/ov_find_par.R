library(ggplot2)
library(tidyverse)
library(patchwork)

# read in male cov and female cov
mcov <- read.table("ov_ERR948801_cov_100k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
mtwon <- median(subset(mcov, chr != "OVOC_OM2")$cov)
mcov <- mcov %>% mutate(mtrans_cov = cov / mtwon)
mcov <- mcov %>% select(!(cov))
mcov <- mcov %>% filter(chr == "OVOC_OM2")

fcov <- read.table("ov_ERR948804_cov_100k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
ftwon <- median(subset(fcov, chr != "OVOC_OM2")$cov)
fcov <- fcov %>% mutate(ftrans_cov = cov / ftwon)
fcov <- fcov %>% select(!(cov))
fcov <- fcov %>% filter(chr == "OVOC_OM2")

# write bed file to find PAR
nonfemale_dominated_bed <- fcov %>% full_join(mcov) %>% mutate(ratio = mtrans_cov/ftrans_cov) %>% filter(ratio > 0.8) %>% select((chr:stop))
write.table(nonfemale_dominated_bed, "ov_nonfemale_dominated.bed", sep="\t", quote = FALSE, col.names=F, row.names=F)
