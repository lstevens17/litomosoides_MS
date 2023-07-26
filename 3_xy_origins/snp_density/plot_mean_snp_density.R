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
di_nigons <- read_tsv("dirofilaria_immitis_location.tsv") %>% filter(query_chr == "dirofilaria_immitis_chrX")

## D IMMITS
# snp
di_SRR13154013_snp <- read.table("di_SRR13154013_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "dirofilaria_immitis_chrX") %>% select(c("chr", "start", "stop", "dens"))
di_SRR13154014_snp <- read.table("di_SRR13154014_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "dirofilaria_immitis_chrX") %>% select(c("chr", "start", "stop", "dens"))
di_SRR13154015_snp <- read.table("di_SRR13154015_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "dirofilaria_immitis_chrX") %>% select(c("chr", "start", "stop", "dens"))

#put all data frames into list
di_list <- list(di_SRR13154013_snp, di_SRR13154014_snp, di_SRR13154014_snp)

#merge all data frames in list
di_df <- di_list %>% reduce(full_join, by=c("chr", "start", "stop")) 
di_df <- di_df %>% mutate(avg=(Reduce("+",.[4:6]))/3)

di_snp_p <- ggplot(data=di_df, aes(x=((start+stop)/2)/1e6, y=avg)) +
  annotate("rect", xmin=13000000/1e6, xmax=28232375/1e6, ymin=0, ymax=10, fill="#d3d3d3", alpha=0.5) +
  geom_line(linewidth=0.5) +
  coord_cartesian(ylim=c(0, 0.02)) +
  theme_bw() +
  facet_wrap(chr ~ .) +
  ggtitle(expression(italic("D. immitis"))) +
  xlab("Position (Mb)") + ylab("Mean male\nheterozygous SNP density") + 
  theme(strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.minor = element_blank()) 

# plot nigons
di_nigon_p <- filter(di_nigons) %>%
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
  theme(legend.position = "None", strip.text = element_blank(), plot.margin = unit(c(0,0.1,0,0.1), "cm"), panel.grid.minor = element_blank()) 

## O VOLVULUS
# snp
ov_ERR676581_snp <- read.table("ov_ERR676581_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))
ov_ERR948794_snp <- read.table("ov_ERR948794_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))
ov_ERR948795_snp <- read.table("ov_ERR948795_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))
ov_ERR948796_snp <- read.table("ov_ERR948796_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))
ov_ERR948800_snp <- read.table("ov_ERR948800_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))
ov_ERR948801_snp <- read.table("ov_ERR948801_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))
ov_ERR948802_snp <- read.table("ov_ERR948802_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "OVOC_OM2") %>% select(c("chr", "start", "stop", "dens"))

#put all data frames into list
ov_list <- list(ov_ERR676581_snp, ov_ERR948794_snp, ov_ERR948795_snp, ov_ERR948796_snp, ov_ERR948800_snp, ov_ERR948801_snp, ov_ERR948802_snp)

#merge all data frames in list
ov_df <- ov_list %>% reduce(full_join, by=c("chr", "start", "stop"), suffix = c("1","2")) 
ov_df <- ov_df %>% mutate(avg=(Reduce("+",.[4:10]))/7)

ov_snp_p <- ggplot(data=ov_df, aes(x=((start+stop)/2)/1e6, y=avg)) +
  annotate("rect", xmin=22100000/1e6, xmax=25485961/1e6, ymin=0, ymax=10, fill="#d3d3d3", alpha=0.5) +
  geom_line(linewidth=0.5) +
  coord_cartesian(ylim=c(0, 0.02)) +
  theme_bw() +
  facet_wrap(chr ~ .) +
  ggtitle(expression(italic("O. volvulus"))) +
  xlab("Position (Mb)") + ylab("Mean male\nheterozygous SNP density") + 
  theme(strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), panel.grid.minor = element_blank()) 

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
  theme(axis.title.y = element_blank(), 
        legend.position = "None", 
        strip.text = element_blank(), 
        plot.margin = unit(c(0,0.1,0,0.1), "cm"), panel.grid.minor = element_blank()) 

## B MALAYI
# snp
bm_SRR3111731_snp <- read.table("bm_SRR3111731_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111738_snp <- read.table("bm_SRR3111738_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111864_snp <- read.table("bm_SRR3111864_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3112012_snp <- read.table("bm_SRR3112012_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111504_snp <- read.table("bm_SRR3111504_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111510_snp <- read.table("bm_SRR3111510_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111544_snp <- read.table("bm_SRR3111544_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111568_snp <- read.table("bm_SRR3111568_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111579_snp <- read.table("bm_SRR3111579_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111581_snp <- read.table("bm_SRR3111581_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111514_snp <- read.table("bm_SRR3111514_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111517_snp <- read.table("bm_SRR3111517_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111629_snp <- read.table("bm_SRR3111629_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111630_snp <- read.table("bm_SRR3111630_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111634_snp <- read.table("bm_SRR3111634_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111636_snp <- read.table("bm_SRR3111636_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111640_snp <- read.table("bm_SRR3111640_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111318_snp <- read.table("bm_SRR3111318_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111319_snp <- read.table("bm_SRR3111319_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR5190290_snp <- read.table("bm_SRR5190290_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR5190289_snp <- read.table("bm_SRR5190289_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111488_snp <- read.table("bm_SRR3111488_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111493_snp <- read.table("bm_SRR3111493_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR3111498_snp <- read.table("bm_SRR3111498_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR5190288_snp <- read.table("bm_SRR5190288_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))
bm_SRR5190287_snp <- read.table("bm_SRR5190287_snp_dens_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens")) %>% filter(chr == "Bm_v4_ChrX_scaffold_001") %>% select(c("chr", "start", "stop", "dens"))

#put all data frames into list
bm_list <- list(bm_SRR3111731_snp,bm_SRR3111738_snp,bm_SRR3111864_snp,bm_SRR3112012_snp,bm_SRR3111504_snp,bm_SRR3111510_snp,bm_SRR3111544_snp,bm_SRR3111568_snp,bm_SRR3111579_snp,bm_SRR3111581_snp,bm_SRR3111514_snp,bm_SRR3111517_snp,bm_SRR3111629_snp,bm_SRR3111630_snp,bm_SRR3111634_snp,bm_SRR3111636_snp,bm_SRR3111640_snp,bm_SRR3111318_snp,bm_SRR3111319_snp,bm_SRR5190290_snp,bm_SRR5190289_snp,bm_SRR3111488_snp,bm_SRR3111493_snp,bm_SRR3111498_snp,bm_SRR5190288_snp,bm_SRR5190287_snp)

#merge all data frames in list
bm_df <- bm_list %>% reduce(full_join, by=c("chr", "start", "stop"), suffix = c("1","2")) 
bm_df <- bm_df %>% mutate(avg=(Reduce("+",.[4:29]))/26) %>% mutate(chr = str_replace(chr, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))

bm_snp_p <- ggplot(data=bm_df, aes(x=((start+stop)/2)/1e6, y=avg)) +
  annotate("rect", xmin=22350000/1e6, xmax=24943668/1e6, ymin=0, ymax=3, fill="#d3d3d3", alpha=0.5) +
  geom_line(linewidth=0.5) +
  coord_cartesian(ylim=c(0, 0.02)) +
  theme_bw() +
  facet_wrap(chr ~ .) +
  ggtitle(expression(italic("B. malayi"))) +
  xlab("Position (Mb)") + ylab("Mean male\nheterozygous SNP density") + 
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

p <- di_snp_p / ov_snp_p / bm_snp_p / di_nigon_p / ov_nigon_p / bm_nigon_p + 
  plot_layout(ncol = 3)

ggsave("mean_snp_density_nigons.png", plot = p, width=10, height=5, units="in")
ggsave("mean_snp_density_nigons.pdf", plot = p, width=10, height=5, units="in")

