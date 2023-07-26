library(ggplot2)
library(tidyverse)
library(patchwork)

# read in locations of nigons
dirofilaria_immitis <- read_tsv("dirofilaria_immitis_location.tsv") 
onchocerca_volvulus <- read_tsv("onchocerca_volvulus_location.tsv")
brugia_malayi <- read_tsv("brugia_malayi_location.tsv") %>% mutate(query_chr = str_replace(query_chr, "Bm_v4_Chr1_scaffold_001", "Bm_v4_Chr1")) %>% 
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
window_size = 1000000

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# plot nigon distrubtion
p1 <- filter(dirofilaria_immitis) %>%
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
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0,0,0), "cm"))  + 
  ggtitle(expression(italic("Dirofilaria immitis")))

p2 <- filter(onchocerca_volvulus) %>%
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
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0,0,0), "cm"))  + 
  ggtitle(expression(italic("Onchocerca volvulus")))

p3 <- filter(brugia_malayi) %>%
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
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0,0,0), "cm"))  + 
  ggtitle(expression(italic("Brugia malayi")))

p4 <- filter(litomosoides_sigmodontis) %>%
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
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0), legend.position = "None", plot.margin = unit(c(0,0,0,0), "cm")) + 
  ggtitle(expression(italic("Litomosoides sigmodontis")))

p <- p1 / p2 / p3 / p4 
ggsave("filarial_nigon_plots.png", plot = p, width=9, height=6, units="in")
ggsave("filarial_nigon_plots.pdf", plot = p, width=9, height=6, units="in")