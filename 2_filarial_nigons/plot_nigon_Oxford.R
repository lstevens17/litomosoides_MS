library(ggplot2)
library(tidyverse)
library(patchwork)
library(forcats)

# nigon set up 
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b", "unassigned" = "grey")

# d_immitis
compare <- read.table("di_vs_ls.oxford_nigon.tsv", col.names=c("buscoID", "nigon", "chr1", "pos1", "chr2", "pos2")) %>% filter(nigon != "unassigned")
compare$chr2 <- factor(compare$chr2, 
                       levels = c("dirofilaria_immitis_chr2", "dirofilaria_immitis_chr3", "dirofilaria_immitis_chr1", "dirofilaria_immitis_chr4", "dirofilaria_immitis_chrX"))
compare$chr1 <- factor(compare$chr1, 
                       levels = c("SUPER_5", "SUPER_X", "SUPER_4", "SUPER_2", "SUPER_3", "SUPER_1"))

di <- compare %>% 
  ggplot(., aes(x=pos2/1e6, y=pos1/1e6, color=nigon)) + 
  geom_point(alpha=0.5, size=2) + facet_grid(chr1~chr2, scales="free", space="free") + theme_bw() + 
  scale_colour_manual(values=cols, name="Nigon\nelement") + xlab(expression(paste("Position in ", italic("D. immitis"), " (Mb)"))) + 
  ylab(expression(paste("Position in ", italic("L. sigmodontis"), " (Mb)")))  +
  theme(panel.spacing=unit(0, "lines"), panel.grid.minor = element_blank()) + 
  scale_x_continuous(expand=c(0,0), breaks=c(5, 10, 15, 20)) + scale_y_continuous(expand=c(0,0), breaks=c(3, 6, 9, 12, 15, 18)) + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.position = c(0.95,0.2), legend.background = element_rect(linetype="solid", colour ="black"))

# o_volvulus
compare <- read.table("ov_vs_ls.oxford_nigon.tsv", col.names=c("buscoID", "nigon", "chr1", "pos1", "chr2", "pos2")) %>% 
  filter(nigon != "unassigned") %>% filter(chr2 == "OVOC_OM1b" | chr2 == "OVOC_OM2" | chr2 == "OVOC_OM3" | chr2 == "OVOC_OM4")
compare$chr2 <- factor(compare$chr2, 
                       levels = c("OVOC_OM1b", "OVOC_OM4", "OVOC_OM3", "OVOC_OM2"))
compare$chr1 <- factor(compare$chr1, 
                       levels = c("SUPER_5", "SUPER_X", "SUPER_2", "SUPER_3", "SUPER_1", "SUPER_4"))

ov <- compare %>% 
  ggplot(., aes(x=pos2/1e6, y=pos1/1e6, color=nigon)) + 
  geom_point(alpha=0.5, size=2) + facet_grid(chr1~chr2, scales="free", space="free") + theme_bw() + 
  scale_colour_manual(values=cols, name="Nigon\nelement") + xlab(expression(paste("Position in ", italic("O. volvulus"), " (Mb)"))) + 
  ylab(expression(paste("Position in ", italic("L. sigmodontis"), " (Mb)")))  +
  theme(panel.spacing=unit(0, "lines"), panel.grid.minor = element_blank()) + 
  scale_x_continuous(expand=c(0,0), breaks=c(5, 10, 15, 20)) + scale_y_continuous(expand=c(0,0), breaks=c(3, 6, 9, 12, 15, 18)) + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.position = c(0.95,0.2), legend.background = element_rect(linetype="solid", colour ="black"))

# b_malayi
compare <- read.table("bm_vs_ls.oxford_nigon.tsv", col.names=c("buscoID", "nigon", "chr1", "pos1", "chr2", "pos2")) %>% 
  filter(nigon != "unassigned") %>% filter(chr2 == "Bm_v4_Chr1_scaffold_001" | chr2 == "Bm_v4_Chr2_contig_001" | chr2 == "Bm_v4_Chr3_scaffold_001" | chr2 == "Bm_v4_Chr4_scaffold_001" | chr2 == "Bm_v4_ChrX_scaffold_001") %>%
  mutate(chr2 = str_replace(chr2, "Bm_v4_Chr1_scaffold_001", "Bm_v4_Chr1")) %>% 
  mutate(chr2 = str_replace(chr2, "Bm_v4_Chr2_contig_001", "Bm_v4_Chr2")) %>% 
  mutate(chr2 = str_replace(chr2, "Bm_v4_Chr3_scaffold_001", "Bm_v4_Chr3")) %>% 
  mutate(chr2 = str_replace(chr2, "Bm_v4_Chr4_scaffold_001", "Bm_v4_Chr4")) %>% 
  mutate(chr2 = str_replace(chr2, "Bm_v4_ChrX_scaffold_001", "Bm_v4_ChrX"))
compare$chr2 <- factor(compare$chr2, 
                       levels = c("Bm_v4_Chr2", "Bm_v4_Chr1", "Bm_v4_Chr3", "Bm_v4_Chr4", "Bm_v4_ChrX"))
compare$chr1 <- factor(compare$chr1, 
                       levels = c("SUPER_4", "SUPER_X", "SUPER_5", "SUPER_1", "SUPER_2", "SUPER_3"))

bm <- compare %>% 
  ggplot(., aes(x=pos2/1e6, y=pos1/1e6, color=nigon)) + 
  geom_point(alpha=0.5, size=2) + facet_grid(chr1~chr2, scales="free", space="free") + theme_bw() + 
  scale_colour_manual(values=cols, name="Nigon\nelement") + xlab(expression(paste("Position in ", italic("B. malayi"), " (Mb)"))) + 
  ylab(expression(paste("Position in ", italic("L. sigmodontis"), " (Mb)")))  +
  theme(panel.spacing=unit(0, "lines"), panel.grid.minor = element_blank()) + 
  scale_x_continuous(expand=c(0,0), breaks=c(5, 10, 15, 20)) + scale_y_continuous(expand=c(0,0), breaks=c(3, 6, 9, 12, 15, 18)) + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=10), legend.position = c(0.95,0.2), legend.background = element_rect(linetype="solid", colour ="black"))

# plot all 
ggsave("di_vs_ls.oxford_nigon.png", plot = di, width = 26, height = 18, units = "cm")
ggsave("ov_vs_ls.oxford_nigon.png", plot = ov, width = 26, height = 18, units = "cm")
ggsave("bm_vs_ls.oxford_nigon.png", plot = bm, width = 26, height = 18, units = "cm")