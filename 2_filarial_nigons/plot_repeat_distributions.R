library(ggplot2)
library(tidyverse)
library(patchwork)

############
## NIGONS ##
############
# read in locations of nigons
dirofilaria_immitis <- read_tsv("dirofilaria_immitis_location.tsv") 
onchocerca_volvulus <- read_tsv("onchocerca_volvulus_location.tsv")
brugia_malayi <- read_tsv("brugia_malayi_location.tsv") %>% 
  mutate(query_chr = str_replace(query_chr, "Bm_v4_Chr1_scaffold_001", "Bm_v4_Chr1")) %>%
  mutate(query_chr = str_replace(query_chr, "Bm_v4_Chr2_contig_001", "Bm_v4_Chr2")) %>%
  mutate(query_chr = str_replace(query_chr, "Bm_v4_Chr3_scaffold_001", "Bm_v4_Chr3")) %>%
  mutate(query_chr = str_replace(query_chr, "Bm_v4_Chr4_scaffold_001", "Bm_v4_Chr4")) %>%
  mutate(query_chr = str_replace(query_chr, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))
litomosoides_sigmodontis <- read_tsv("litomosoides_sigmodontis_location.tsv")

# set order 
dirofilaria_immitis$query_chr <- 
  factor(dirofilaria_immitis$query_chr, levels=c("dirofilaria_immitis_chr2", "dirofilaria_immitis_chr3", "dirofilaria_immitis_chr1", "dirofilaria_immitis_chr4", "dirofilaria_immitis_chrX"))
onchocerca_volvulus$query_chr <- 
  factor(onchocerca_volvulus$query_chr, levels=c("OVOC_OM1b", "OVOC_OM4", "OVOC_OM3", "OVOC_OM2"))
brugia_malayi$query_chr <- 
  factor(brugia_malayi$query_chr, levels=c("Bm_v4_Chr3", "Bm_v4_Chr2", "Bm_v4_Chr1", "Bm_v4_Chr4", "Bm_v4_ChrX")) 
litomosoides_sigmodontis$query_chr <- 
  factor(litomosoides_sigmodontis$query_chr, levels=c("SUPER_1", "SUPER_3", "SUPER_2", "SUPER_5", "SUPER_4", "SUPER_X"))

# set up nigon window size and colour
window_size = 500000

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# plot nigon distrubtion
di_nigon_p <- filter(dirofilaria_immitis) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(. ~ query_chr, scales = "free_x", space = "free_x") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO\ncount (n)") + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "None") 

ov_nigon_p <- filter(onchocerca_volvulus) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 200) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(. ~ query_chr, scales = "free_x", space = "free_x") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO\ncount (n)") + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "None") 

bm_nigon_p <- filter(brugia_malayi) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(. ~ query_chr, scales = "free_x", space = "free_x") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO\ncount (n)") + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        panel.grid.minor = element_blank(),
        legend.position = "None") 

ls_nigon_p <- filter(litomosoides_sigmodontis) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(. ~ query_chr, scales = "free_x", space = "free_x") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO\ncount (n)")  + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "None") 


#############
## REPEATS ##
#############
# read in repeats
ls_repeats <- read.table("nxLitSigm11.1.primary_red_200kb.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr != "scaffold_MT_1")
di_repeats <- read.table("dimmitis_WSI_2.2_red_200kb.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens"))
ov_repeats <- read.table("onchocerca_volvulus.PRJEB513.WBPS18.genomic_red_200kb.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens"))
bm_repeats <- read.table("brugia_malayi.PRJNA10729.WBPS18.genomic_red_200kb.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>%
  mutate(chr = str_replace(chr, "Bm_v4_Chr1_scaffold_001", "Bm_v4_Chr1")) %>%
  mutate(chr = str_replace(chr, "Bm_v4_Chr2_contig_001", "Bm_v4_Chr2")) %>%
  mutate(chr = str_replace(chr, "Bm_v4_Chr3_scaffold_001", "Bm_v4_Chr3")) %>%
  mutate(chr = str_replace(chr, "Bm_v4_Chr4_scaffold_001", "Bm_v4_Chr4")) %>%
  mutate(chr = str_replace(chr, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))

# set order 
di_repeats$chr <- 
  factor(di_repeats$chr, levels=c("dirofilaria_immitis_chr2", "dirofilaria_immitis_chr3", "dirofilaria_immitis_chr1", "dirofilaria_immitis_chr4", "dirofilaria_immitis_chrX"))
ov_repeats$chr <- 
  factor(ov_repeats$chr, levels=c("OVOC_OM1b", "OVOC_OM4", "OVOC_OM3", "OVOC_OM2"))
bm_repeats$chr <- 
  factor(bm_repeats$chr, levels=c("Bm_v4_Chr3", "Bm_v4_Chr2", "Bm_v4_Chr1", "Bm_v4_Chr4", "Bm_v4_ChrX")) 
ls_repeats$chr <- 
  factor(ls_repeats$chr, levels=c("SUPER_1", "SUPER_3", "SUPER_2", "SUPER_5", "SUPER_4", "SUPER_X"))

# make plots
ls_rep_p <- ggplot(ls_repeats, aes(x=((start+stop)/2)/1e6, y=dens)) + 
  geom_point(alpha=0.2) + 
  geom_smooth(span = 0.2, color="black", se=F) +
  coord_cartesian(ylim=c(0,1)) + 
  facet_grid(.~chr, space="free_x", scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Repeat density")  + 
  ggtitle(expression(italic("L. sigmodontis"))) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.minor = element_blank()) 

di_rep_p <- ggplot(di_repeats, aes(x=((start+stop)/2)/1e6, y=dens)) + 
  geom_point(alpha=0.2) + 
  geom_smooth(span = 0.2, color="black", se=F) +
  coord_cartesian(ylim=c(0,1)) + 
  facet_grid(.~chr, space="free_x", scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Repeat density") + 
  ggtitle(expression(italic("D. immitis"))) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.minor = element_blank()) 

ov_rep_p <- ggplot(ov_repeats, aes(x=((start+stop)/2)/1e6, y=dens)) + 
  geom_point(alpha=0.2) + 
  geom_smooth(span = 0.2, color="black", se=F) +
  coord_cartesian(ylim=c(0,1)) + 
  facet_grid(.~chr, space="free_x", scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Repeat density") + 
  ggtitle(expression(italic("O. volvulus"))) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.minor = element_blank()) 

bm_rep_p <- ggplot(bm_repeats, aes(x=((start+stop)/2)/1e6, y=dens)) + 
  geom_point(alpha=0.2) + 
  geom_smooth(span = 0.2, color="black", se=F) +
  coord_cartesian(ylim=c(0,1)) + 
  facet_grid(.~chr, space="free_x", scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Repeat density") + 
  ggtitle(expression(italic("B. malayi"))) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.minor = element_blank()) 

p <- ls_rep_p / ls_nigon_p / di_rep_p / di_nigon_p / ov_rep_p / ov_nigon_p / bm_rep_p / bm_nigon_p

ggsave("repeat_distribution.png", plot = p, width=10, height=12, units="in")
ggsave("repeat_distribution.pdf", plot = p, width=10, height=12, units="in")
