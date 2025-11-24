#Filter Depth and MAF QC
# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# Read the imiss file from vcftools
imiss <- read_tsv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/New_Filter_QC/missing_ind.imiss", comment = "#")

# Basic summary
summary_stats <- imiss %>%
  summarize(n_indiv = n(),
            mean_fmiss = mean(F_MISS, na.rm=TRUE),
            median_fmiss = median(F_MISS, na.rm=TRUE),
            p95_fmiss = quantile(F_MISS, 0.95, na.rm=TRUE))
print(summary_stats)

# Flag outliers, e.g., F_MISS > 0.10 (10%)
threshold <- 0.10
flagged <- imiss %>% filter(F_MISS > threshold)
write.table(flagged$INDV, file="high_missing_individuals.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# Histogram of per-individual missingness
ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
  labs(x = "Fraction missing per individual (F_MISS)", y = "Count",
       title = "Per-individual missingness (vcftools .imiss)")

# Density plot
ggplot(imiss, aes(x = F_MISS)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
  labs(x = "Fraction missing per individual (F_MISS)", y = "Density",
       title = "Density of per-individual missingness")

# Bar plot of worst individuals (top 20 by F_MISS)
topN <- imiss %>% arrange(desc(F_MISS)) %>% slice_head(n = 20)
ggplot(topN, aes(x = reorder(INDV, F_MISS), y = F_MISS)) +
  geom_col(fill = "tomato") +
  coord_flip() +
  labs(x = "Individual", y = "F_MISS",
       title = "Top 20 individuals by missingness")
