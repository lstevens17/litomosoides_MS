library(ggplot2)
library(tidyverse)
library(patchwork)

telomere_bed <- read.table("nxLitSigm11.1.primary.fa_telomere_density_1kb.bed", col.names=c("chr", "start", "end", "count", "span", "window", "dens")) %>% filter(chr != "scaffold_MT_1")

p <- ggplot(data=telomere_bed, aes(x=start/1e6, y=count)) + 
  geom_line() + 
  facet_wrap(~chr, scales="free_x", nrow=2) + 
  theme_bw() + xlab("Position (Mb)") + ylab("Counts of TTAGGC\nper 1kb window") 

ggsave("telomere_distributions.png", plot = p, width=10, height=5, units="in")
ggsave("telomere_distributions.pdf", plot = p, width=10, height=5, units="in")
