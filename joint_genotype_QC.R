library(tidyverse)

# Load heterozygosity
het <- read_table("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/het_merged_qc.het", col_types = cols())

het <- het |>
  rename(IID = `#IID`) |>
  mutate(
    obs_het = 1 - (`O(HOM)` / OBS_CT),   # observed heterozygosity
    exp_het = 1 - (`E(HOM)` / OBS_CT)    # expected heterozygosity
  ) |> 
  filter(!IID %in% c("GAL-057","GAL-127","GAL-225")) #remove outliers

### 2. Read population metadata and join
# pop_map.txt: at least two columns: IID, population
pop_map <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_no_57_127_pca_env_data.csv") |> 
  filter(!IID %in% c("GAL-057","GAL-127","GAL-225"))


het_pop <- het |>
  left_join(pop_map, by = "IID")

### 3. Boxplot of observed heterozygosity by population
ggplot(het_pop, aes(x = Cluster, y = obs_het, fill = Cluster)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Population",
    y = "Observed heterozygosity",
    title = "Per-sample observed heterozygosity by population"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### 4. Optional: overlay points on boxplot
ggplot(het_pop, aes(x = population, y = obs_het, fill = population)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  theme_bw() +
  labs(
    x = "Population",
    y = "Observed heterozygosity",
    title = "Per-sample observed heterozygosity by population"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### 5. Optional: plot inbreeding coefficient F by population
ggplot(het_pop, aes(x = population, y = F, fill = population)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Population",
    y = "Inbreeding coefficient (F)",
    title = "Per-sample F by population"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


library(tidyverse)

miss <- read_table("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/final_merged_qc.smiss", col_types = cols()) |>
  rename(IID = `#IID`)

het_miss <- het_pop |>
  left_join(miss, by = "IID")

ggplot(het_miss, aes(x = F_MISS, y = obs_het, color = Cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_bw() +
  labs(
    x = "Per-sample missingness",
    y = "Observed heterozygosity",
    color = "Population",
    title = "Heterozygosity vs missingness by population"
  )


#Filtered VCF file No samples 57 and 127
#read eigenvectors and eigenvalues (Produced using PLINK2)

eigenValues <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/combined_gal_mainland_snp_pca_results.eigenval", delim = " ", col_names = F)

eigenVectors <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/combined_gal_mainland_snp_pca_results.eigenvec", delim = "\t", col_names = T) |>
  rename(IID = `#IID`)

#env_data <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/sample_population_info.csv")

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

eigen_env <- left_join(het_pop, eigenVectors, by = "IID")

#Plot PCA
ggplot(data = eigen_env,mapping = aes(x = PC1, y = PC2, color = Bioregion, shape = Cluster)) +
  geom_point(size = 3, show.legend = T ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Gal & MAinland - No Outliers (57 & 127)",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()


#############################
############################
######Het per VCF###########
library(tidyverse)

### Read heterozygosity files
het_gal <- read_table("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/gal_het_original_maxDP_60.het", col_types = cols()) |>
  rename(IID = `#IID`) |>
  mutate(
    obs_het = 1 - (`O(HOM)` / OBS_CT),
    exp_het = 1 - (`E(HOM)` / OBS_CT),
    batch = "Galapagos"
  ) |> 
  filter(!IID %in% c("GAL-057","GAL-127","GAL-131"))

het_mainland <- read_table("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/mainland_het_original_maxDP_95.het", col_types = cols()) |>
  rename(IID = `#IID`) |>
  mutate(
    obs_het = 1 - (`O(HOM)` / OBS_CT),
    exp_het = 1 - (`E(HOM)` / OBS_CT),
    batch = "Mainland"
  )|> 
  filter(!IID %in% c("GAL-225"))

### Combine
het_combined <- bind_rows(het_gal, het_mainland)

### Read population metadata
#pop_map <- read_tsv("population_map.txt", col_types = cols())
# Expected columns: IID, population

het_pop <- het_combined |>
  left_join(pop_map, by = "IID")

### Read missingness (optional QC check)#####
miss_gal <- read_table("gal_miss_original.smiss", col_types = cols()) |>
  rename(FID = `#FID`) |>
  mutate(batch = "Galapagos")

miss_mainland <- read_table("mainland_miss_original.smiss", col_types = cols()) |>
  rename(FID = `#FID`) |>
  mutate(batch = "Mainland")

miss_combined <- bind_rows(miss_gal, miss_mainland)

het_pop <- het_pop |>
  left_join(miss_combined |> select(IID, F_MISS), by = "IID")

### Boxplot by population################
ggplot(het_pop, aes(x = Cluster, y = obs_het, fill = Cluster)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.7) +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  theme_bw(base_size = 14) +
  labs(
    x = "Population",
    y = "Observed heterozygosity",
    title = "DP filtered Per-sample heterozygosity"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### Het vs missingness (to check if artifact persists)######
ggplot(het_pop, aes(x = F_MISS, y = obs_het, color = population)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    x = "Per-sample missingness (within batch)",
    y = "Observed heterozygosity",
    color = "Population",
    title = "Heterozygosity vs missingness (within-batch calculation)"
  )

### Summary statistics by population
het_pop |>
  group_by(Cluster) |>
  summarise(
    n = n(),
    mean_het = mean(obs_het, na.rm = TRUE),
    sd_het = sd(obs_het, na.rm = TRUE),
    median_het = median(obs_het, na.rm = TRUE),
    mean_miss = mean(F_MISS, na.rm = TRUE),
    mean_F = mean(F, na.rm = TRUE),
    .groups = "drop"
  )
#############################
############################

# From your het_pop table
ggplot(het_pop, aes(x = Cluster, y = F, fill = Cluster)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Inbreeding coefficient by population", y = "F")
#############################




######Depth per VCF###########
library(tidyverse)

gal_dp <- read_tsv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/gal_mean_dp.tsv", col_names = c("IID", "mean_DP")) |>
  mutate(batch = "Galapagos")

mainland_dp <- read_tsv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/mainland_mean_dp.tsv", col_names = c("IID", "mean_DP")) |>
  mutate(batch = "Mainland")

dp_all <- bind_rows(gal_dp, mainland_dp)

dp_all %>%
  group_by(batch) %>%
  summarise(
    n = n(),
    mean_dp = mean(mean_DP),
    sd_DP   = sd(mean_DP),
    min_DP  = min(mean_DP),
    max_DP  = max(mean_DP)
  )



