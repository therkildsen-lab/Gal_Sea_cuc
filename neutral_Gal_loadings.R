suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(grid)
  library(tidyverse)
})

# Read SNP loadings for PC1
snp_loadings <- fread("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/Neutral_Gal_Only/full_neutral_genome_snp_plink_biallelic_snp_pca.eigenvec.var", header = TRUE) |>
  select(`#CHROM`,ID,PC1) |>
  mutate(PC1_sqrd = PC1^2)
# #
# # #prepare data
data <- snp_loadings %>%
  filter(is.finite(PC1_sqrd)) %>%
  mutate(Rank_unique = dense_rank(desc(PC1_sqrd))) %>%  # 1 = highest value; ties share rank, many loadings are similar
  mutate(
    Top_100 = Rank_unique <= 100,
    Top_0_1_Percent = row_number(desc(PC1_sqrd)) <= ceiling(0.001 * n())
  ) %>%
  arrange(Rank_unique, desc(PC1_sqrd)) %>%
  filter(Top_0_1_Percent) %>% # Filter for top-ranked SNPs (top 0.1%)
  mutate(POS = sub("^[^:]*:([^:]+).*", "\\1", ID)) # position from your ID field
# #
# # #save clean DF
write_csv(data,"Data/neutral_top_snps_0_1.csv")

#save top_SNPs IDs for furhter PLINK analyses
top_snps <- unique(data$ID)
# 
# # Write SNP list for PLINK2 extraction
write.table(top_snps, file="neutral_top0.1pct_snps.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

#load data
#top_snps_1 <- read_csv("Data/top_snps_0_1.csv",show_col_types = FALSE)

# 1) Normalize minimal types, compute offsets once, cumulative position, and axis centers
genome_df <- data %>%
  transmute(
    CHR_NUM = as.integer(`#CHROM`),           # keep numbers as provided (non-human allowed)
    POS     = as.numeric(POS),                  
    PC1,
    PC1_sqrd,
    Top_100 = as.logical(Top_100)
  ) %>%
  filter(!is.na(CHR_NUM), !is.na(POS)) %>%
  arrange(CHR_NUM, POS)

chr_info <- genome_df %>%
  group_by(CHR_NUM) %>%
  summarise(chr_len = max(POS, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR_NUM) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

df_cum <- genome_df %>%
  left_join(chr_info, by = "CHR_NUM") %>%
  mutate(
    pos_cum  = POS + offset,
    bg_group = CHR_NUM %% 2           # alternate two colors by chromosome
  )

axis_df <- df_cum %>%
  group_by(CHR_NUM) %>%
  summarise(center = (min(pos_cum) + max(pos_cum)) / 2, .groups = "drop")

# 2) Plot PC1 along genome (alternating grey/blue)
p_pc1 <- ggplot(df_cum[df_cum$PC1>0,], aes(x = pos_cum, y = PC1)) +
  geom_point(aes(color = factor(bg_group)), alpha = 0.6, size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = c("0" = "grey70", "1" = "#2C7FB8")) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR_NUM,
                     expand = expansion(mult = c(0.001, 0.01))) +
  labs(title = "Genome-wide Neutral SNP subset of 0.1% loadings for PC1 (N=23635)",
       x = "Chromosomes", y = "PC1 loading") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

p_pc1
