library(ggplot2)
library(tidyverse)
library(patchwork)

####################
### NIGON SET UP ###
####################
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")
window_size = 500000

# parse nigon files files 
ov_nigons <- read_tsv("onchocerca_volvulus_location.tsv") %>% filter(query_chr == "OVOC_OM2")
bm_nigons <- read_tsv("brugia_malayi_location.tsv") %>% filter(query_chr == "Bm_v4_ChrX_scaffold_001")

ov_y_contigs <- read.table("onchocerca_volvulus_y_contigs_vs_genome.filtered_50kb_windows.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))

ov_y_p <- ggplot(data=ov_y_contigs, aes(x=((start+stop)/2)/1e6, y=dens)) +
  annotate("rect", xmin=22100000/1e6, xmax=25485961/1e6, ymin=0, ymax=2, fill="#d3d3d3", alpha=1) +
  geom_line(linewidth=0.5) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw() +
  facet_wrap(chr ~ .) +
  ggtitle(expression(italic("O. volvulus"))) +
  xlab("Position (Mb)") + ylab("Y contig alignment coverage") + 
  theme(strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.minor = element_blank()) 

# plot nigons
ov_nigon_p <- filter(ov_nigons) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-(window_size/2))/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_wrap(query_chr ~ ., nrow=1, , scales="free_x") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") + 
  theme(legend.position = "None", 
        strip.text = element_blank(), 
        plot.margin = unit(c(0,0.1,0,0.1), "cm"), panel.grid.minor = element_blank()) 



bm_y_contigs <- read.table("brugia_malayi_y_contigs_vs_genome.filtered_50kb_windows.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens")) %>%
  mutate(chr = str_replace(chr, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))

bm_y_p <- ggplot(data=bm_y_contigs, aes(x=((start+stop)/2)/1e6, y=dens)) +
  annotate("rect", xmin=22350000/1e6, xmax=24943668/1e6, ymin=0, ymax=2, fill="#d3d3d3", alpha=1) +
  geom_line(linewidth=0.5) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw() +
  facet_wrap(chr ~ .) +
  ggtitle(expression(italic("B. malayi"))) +
  xlab("Position (Mb)") + 
  theme(strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), panel.grid.minor = element_blank()) 

# plot nigons
bm_nigon_p <- filter(bm_nigons) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-(window_size/2))/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_wrap(query_chr ~ ., nrow=1, , scales="free_x") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") + 
  theme(axis.title.y = element_blank(), legend.position = "None", strip.text = element_blank(), plot.margin = unit(c(0,0.1,0,0.1), "cm"), panel.grid.minor = element_blank()) 

p <- ov_y_p / bm_y_p / ov_nigon_p / bm_nigon_p + 
  plot_layout(ncol = 2)

ggsave("y_contig_coverage.png", plot = p, width=9, height=6, units="in")
ggsave("y_contig_coverage.pdf", plot = p, width=9, height=6, units="in")
