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
di_nigons <- read_tsv("dirofilaria_immitis_location.tsv") %>% filter(query_chr == "dirofilaria_immitis_chrX")
ov_nigons <- read_tsv("onchocerca_volvulus_location.tsv") %>% filter(query_chr == "OVOC_OM2")
bm_nigons <- read_tsv("brugia_malayi_location.tsv") %>% filter(query_chr == "Bm_v4_ChrX_scaffold_001")  %>%
  mutate(query_chr = str_replace(query_chr, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))

#############
# DIMMITIS ##
#############

# read in male cov and female cov
mcov <- read.table("di_SRR13154015_cov_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
mtwon <- median(subset(mcov, chr != "dirofilaria_immitis_chrX")$cov)
mcov <- mcov %>% mutate(trans_cov = cov / mtwon)
mcov <- mcov %>% filter(chr == "dirofilaria_immitis_chrX")

fcov <- read.table("di_ERR034941_cov_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
ftwon <- median(subset(fcov, chr != "dirofilaria_immitis_chrX")$cov)
fcov <- fcov %>% mutate(trans_cov = cov / ftwon)
fcov <- fcov %>% filter(chr == "dirofilaria_immitis_chrX")

mcov <- mcov %>% mutate(sex = "male")
fcov <- fcov %>% mutate(sex = "female")

bothcov <- rbind(mcov, fcov)

# plot cov
di_cov_p <- ggplot(bothcov) + 
  annotate("rect", xmin=13000000/1e6, xmax=28232375/1e6, ymin=0, ymax=3, fill="#d3d3d3", alpha=1) +
  geom_line(aes(x=((start+stop)/2)/1e6, y=trans_cov*2, colour=sex), linewidth=0.9, alpha=1) +
  theme_bw() + 
  facet_wrap(chr ~ . ) + 
  xlab("Position (Mb)") + ylab("Normalised coverage (N)") +
  ylim(0, 3) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "None",
        plot.title = element_text(face = "bold")
  ) +
  ggtitle(expression(italic("Dirofilaria immitis"))) +
  scale_colour_manual(values=c("#0092c1", "black"))

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
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(legend.position = "None", strip.text = element_blank())

##############
# OVOLVULUS ##
##############

# read in male cov and female cov
mcov <- read.table("ov_ERR948801_cov_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
mtwon <- median(subset(mcov, chr != "OVOC_OM2")$cov)
mcov <- mcov %>% mutate(trans_cov = cov / mtwon)
mcov <- mcov %>% filter(chr == "OVOC_OM2")

fcov <- read.table("ov_ERR948804_cov_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
ftwon <- median(subset(fcov, chr != "OVOC_OM2")$cov)
fcov <- fcov %>% mutate(trans_cov = cov / ftwon)
fcov <- fcov %>% filter(chr == "OVOC_OM2")

mcov <- mcov %>% mutate(sex = "male")
fcov <- fcov %>% mutate(sex = "female")

bothcov <- rbind(mcov, fcov)

# plot cov
ov_cov_p <- ggplot(bothcov) + 
  annotate("rect", xmin=22100000/1e6, xmax=25485961/1e6, ymin=0, ymax=3, fill="#d3d3d3", alpha=1) +
  geom_line(aes(x=((start+stop)/2)/1e6, y=trans_cov*2, colour=sex), linewidth=0.9, alpha=1) +
  theme_bw() + 
  facet_wrap(chr ~ . ) + 
  xlab("Position (Mb)") + ylab("Normalised coverage (N)") +
  ylim(0, 3) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "None",
        plot.title = element_text(face = "bold")
        ) +
  ggtitle(expression(italic("Onchocerca volvulus"))) +
  scale_colour_manual(values=c("#0092c1", "black"))

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
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(legend.position = "None", strip.text = element_blank())

############
# BMALAYI ##
############

# read in male cov and female cov
mcov <- read.table("bm_SRR3111510_cov_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
mtwon <- median(subset(mcov, chr != "Bm_v4_ChrX_scaffold_001")$cov)
mcov <- mcov %>% mutate(trans_cov = cov / mtwon)
mcov <- mcov %>% filter(chr == "Bm_v4_ChrX_scaffold_001")  %>%
  mutate(chr = str_replace(chr, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))

fcov <- read.table("bm_SRR488256_cov_50k_windows.regions.bed", col.names=c("chr", "start", "stop", "cov")) 
ftwon <- median(subset(fcov, chr != "Bm_v4_ChrX_scaffold_001")$cov)
fcov <- fcov %>% mutate(trans_cov = cov / ftwon)
fcov <- fcov %>% filter(chr == "Bm_v4_ChrX_scaffold_001")  %>%
  mutate(chr = str_replace(chr, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))

mcov <- mcov %>% mutate(sex = "male")
fcov <- fcov %>% mutate(sex = "female")

bothcov <- rbind(mcov, fcov)

# plot cov
bm_cov_p <- ggplot(bothcov) + 
  annotate("rect", xmin=22350000/1e6, xmax=24943668/1e6, ymin=0, ymax=3, fill="#d3d3d3", alpha=1) +
  geom_line(aes(x=((start+stop)/2)/1e6, y=trans_cov*2, colour=sex), linewidth=0.9, alpha=1) +
  theme_bw() + 
  facet_wrap(chr ~ . ) + 
  xlab("Position (Mb)") + ylab("Normalised coverage (N)") +
  ylim(0, 3) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "None", 
        plot.title = element_text(face = "bold")
        ) + 
  ggtitle(expression(italic("Brugia malayi"))) +
  scale_colour_manual(values=c("#0092c1", "black"))

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
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(legend.position = "None", strip.text = element_blank())

############
### PLOT ###
############
# assemble
p <- di_cov_p / ov_cov_p / bm_cov_p / di_nigon_p / ov_nigon_p / bm_nigon_p +
  plot_layout(nrow=2, ncol=3)

ggsave("Figure3_PAR_coverage_nigons.png", plot = p, width=12, height=6, units="in")
ggsave("Figure3_PAR_coverage_nigons.pdf", plot = p, width=12, height=6, units="in")


