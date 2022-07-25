setwd("D:/Users/Chauncey/Desktop/16S/SZW20220406/Analysis/7_OTUs_alpha/Total")

library(tidyverse)
library(ggsignif)
library(patchwork)

sample.list <- read_tsv("sample list.txt",
                        col_names = FALSE)
data <- read_tsv("alpha-summary.tsv")
colnames(data)[1] <- "SampleID"
colnames(data)[2] <- "shannon"
colnames(data)[5] <- "chao1"
# data <- read_tsv("Total.final.otu.chao1.txt")
# colnames(data)[1] <- "SampleID"

group <- read_tsv("../../Total-color3.txt") %>% 
  filter(SampleID %in% sample.list$X1)

# write_tsv(group, "../../Total-color4.txt")

mydata <- right_join(data, group) %>% 
  mutate(
    Group = factor(
      Group,
      levels = c("Control", "L.fermentum", "L.animalis", "L.salvarius", "L.reuteri")
    )
  ) %>% 
  arrange(Group)

mycol <- unique(mydata$Color)

myCompare <- list(
  c("Control", "L.fermentum"),
  c("Control", "L.animalis"),
  c("Control", "L.salvarius"),
  c("Control", "L.reuteri")
)

p_sd <- ggplot(
  mydata,
  aes(Group, chao1)
) +
  geom_point(
    data = mydata,
    aes(Group, chao1, color = Group),
    size = 2.4,
    position = position_jitterdodge(0.6),
    alpha = .3
  ) +
  geom_signif(
    comparisons = myCompare,
    test = "wilcox.test",
    step_increase = .18,
    textsize = 4.5,
    map_signif_level = TRUE
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    size = 1.2,
    width = .35,
    aes(color = Group)
  ) +
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 5,
    aes(color = Group)
  ) +
  scale_color_manual(values = mycol) +
  ggtitle("") +
  ylab("Chao index") +
  # scale_x_discrete(labels = c("T0", "T2", "T3", "T0/T3")) +
  scale_y_continuous(limits = c(7000, 15500)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(colour = "black", size = 18),
    legend.position = "none",
    plot.title = element_text(hjust = .5, size = 26),
    axis.title.y = element_text(size = 24, margin = unit(c(0, 0.5, 0, 0.5), "lines")),
    axis.text.x = element_text(size = 18, angle = 30, hjust = 1, face = 3),
    panel.background = element_rect(fill = "white", color = "white", size = 4),
    panel.border = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(filename = "chao_sd_0725.pdf", plot = p_sd, height = 4, width = 6)



# 2 -----------------------------------------------------------------------

# sample.list <- read_tsv("sample list.txt",
#                         col_names = FALSE)
# data <- read_tsv("Total.final.otu.shannon.txt")
# colnames(data)[1] <- "SampleID"
# 
# group <- read_tsv("../../Total-color2.txt") %>% 
#   filter(SampleID %in% sample.list$X1)
# 
# # write_tsv(group, "../../Total-color3.txt")
# 
# mydata <- right_join(data, group) %>% 
#   mutate(
#     Group = factor(
#       Group,
#       levels = c("Control", "L.fermentum", "L.animalis", "L.salvarius", "L.reuteri")
#     )
#   ) %>% 
#   arrange(Group)

# mycol <- unique(mydata$Color)

# myCompare <- list(
#   c("Control", "L.fermentum"),
#   c("Control", "L.animalis"),
#   c("Control", "L.salvarius"),
#   c("Control", "L.reuteri")
# )

p_sd2 <- ggplot(
  mydata,
  aes(Group, shannon)
) +
  geom_point(
    data = mydata,
    aes(Group, shannon, color = Group),
    size = 2.4,
    position = position_jitterdodge(0.6),
    alpha = .3
  ) +
  geom_signif(
    comparisons = myCompare,
    test = "wilcox.test",
    step_increase = .18,
    textsize = 4.5,
    map_signif_level = TRUE
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    size = 1.2,
    width = .35,
    aes(color = Group)
  ) +
  stat_summary(
    fun.y = mean,
    geom = "point",
    size = 5,
    aes(color = Group)
  ) +
  scale_color_manual(values = mycol) +
  ggtitle("") +
  ylab("Shannon index") +
  # scale_x_discrete(labels = c("T0", "T2", "T3", "T0/T3")) +
  scale_y_continuous(limits = c(5, 7.5),
                     breaks = c(5, 6, 7)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(colour = "black", size = 18),
    legend.position = "none",
    plot.title = element_text(hjust = .5, size = 26),
    axis.title.y = element_text(size = 24, margin = unit(c(0, 0.5, 0, 0.5), "lines")),
    axis.text.x = element_text(size = 18, angle = 30, hjust = 1, face = 3),
    panel.background = element_rect(fill = "white", color = "white", size = 4),
    panel.border = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(filename = "shannon_sd_0725.pdf", plot = p_sd2, height = 4, width = 6)

pp <- p_sd + p_sd2
ggsave(filename = "combined_sd_0725.pdf", plot = pp, height = 4, width = 8)
