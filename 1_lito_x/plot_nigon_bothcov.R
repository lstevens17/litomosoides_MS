library(ggplot2)
library(tidyverse)
library(patchwork)
library(png)

# read in Hi-C image
hic_image <- readPNG("2023.07.19.13.50.40.HiCImage.png", native = TRUE) 

# read in cov bed (15 male, 16 female)
nxLitSigm15 <- read.table("nxLitSigm15_hic_10perc_100k_windows.regions.bed", 
                          col.names=c("chr", "start", "stop", "cov")) %>% filter(chr != "scaffold_MT_1")
nxLitSigm16 <- read.table("nxLitSigm16_hic_10perc_100k_windows.regions.bed", 
                          col.names=c("chr", "start", "stop", "cov")) %>% filter(chr != "scaffold_MT_1")

# normalise coverage 
mtwon <- median(subset(nxLitSigm15, chr != "SUPER_X")$cov)
mcov <- nxLitSigm15 %>% mutate(trans_cov = cov / mtwon)

ftwon <- median(subset(nxLitSigm16, chr != "SUPER_X")$cov)
fcov <- nxLitSigm16 %>% mutate(trans_cov = cov / ftwon)

mcov <- mcov %>% mutate(sex = "male")
fcov <- fcov %>% mutate(sex = "female")

bothcov <- rbind(mcov, fcov)

# read in locations of nigons
nigon_location <- read_tsv("nxLitSigm11_location.tsv")

# set up nigon window size and colour
window_size = 1000000

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# plot nigon distrubtion
p1 <- filter(nigon_location) %>%
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
  facet_wrap(query_chr ~ ., ncol=6) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() + 
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(strip.text.y = element_text(angle = 0), 
        legend.position = "None", 
        plot.margin = unit(c(0,0,0,0), "cm"), 
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        ) 

# plot read cov distrubtion
p2 <- ggplot(bothcov) + 
  geom_line(aes(x=((start+stop)/2)/1e6, y=trans_cov*2, colour=sex), linewidth=0.9, alpha=0.75) +
  facet_wrap(~chr, ncol=6) + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab(expression(paste("Normalised coverage (", italic("N"), ")"))) + 
  scale_colour_manual(values=c("#0092c1", "black"), name="Sex") + 
  ylim(0, 2.5) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.06, 0.22), 
        legend.background = element_rect(linetype="solid", colour ="black", linewidth=0.25),
        legend.key.height= unit(0.2, 'cm'), 
        legend.key.width= unit(0.2, 'cm'))

# assemble plot
p <- wrap_elements(hic_image) + (p2 / p1) + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 18, face="bold"))

# save
ggsave("Figure1_hic_nigon_bothcov.pdf", plot = p, width=15, height=5, units="in")
ggsave("Figure1_hic_nigon_bothcov.png", plot = p, width=15, height=5, units="in")

