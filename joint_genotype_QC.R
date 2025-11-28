library(tidyverse)

### 2. Read population metadata and join
# pop_map.txt: at least two columns: IID, population
pop_map <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_pop_data.csv") 

#MISSINGNESS

miss <- read_table("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/final_joint_qc.smiss", col_types = cols()) 
het_miss <- pop_map |>
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

# Load heterozygosity
het <- read_table("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/het_joint_qc.het", col_types = cols())

het <- het |>
  
  mutate(
    obs_het = 1 - (`O(HOM)` / OBS_CT),   # observed heterozygosity
    exp_het = 1 - (`E(HOM)` / OBS_CT)    # expected heterozygosity
  ) |> 
  #filter(!IID %in% c("GAL-057","GAL-127","GAL-225")) #remove outliers
filter(!IID %in% c("GAL-225")) #remove duplicate sample



het_pop <- het |>
  left_join(pop_map, by = "IID") |> 
  mutate(
    new_Cluster = factor(new_Cluster,
                     levels = c("Mainland","LEFT", "Center", "RIGHT",  "Outlier")))



### 4. Optional: overlay points on boxplot
ggplot(het_pop, aes(x = new_Cluster, y = obs_het, fill = new_Cluster)) +
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
ggplot(het_pop, aes(x = new_Cluster, y = F, fill = new_Cluster)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
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


#No outlier het and inbreeding (however not recalculating SNPs in PLINK because then it is not comparable) So just filtering here in R
het_no_outlier <- het_pop |> 
  filter(!IID %in% c("GAL-057","GAL-127")) |> 
mutate(
    new_Cluster = factor(new_Cluster,
                     levels = c("Mainland","LEFT", "Center", "RIGHT"))
    # or:
    # Cluster = fct_relevel(Cluster, "LEFT", "Center", "RIGHT", "Outlier")
  )

### 3. Boxplot of observed heterozygosity by population

my_palette1 <- c(
  Mainland = "#FF61C3",
  LEFT = "#CFA127",
  Center = "#56B1F7",
  RIGHT = "#62D99C"
)

my_palette2 <- c(
  RIGHT = "#53051D",
  Mainland = "#377EB8",
  Center = "#F12761",
  LEFT = "#F2CE1B"
)

ggplot(het_no_outlier, aes(x = new_Cluster, y = obs_het, fill = new_Cluster)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = my_palette2) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  theme_bw() +
  labs( 
    x = "Population",
    y = "Observed heterozygosity",
    title = "Per-sample observed heterozygosity by PCA cluster"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### 5. Optional: plot inbreeding coefficient F by population with ordered populations: LEFT, Center, RIGHT

ggplot(het_no_outlier, aes(x = new_Cluster, y = F, fill = new_Cluster)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = my_palette2) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Population",
    y = "Inbreeding coefficient (F)",
    title = "Per-sample F by PCA cluster"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

    



















#PCA
#read eigenvectors and eigenvalues (Produced using PLINK2)

eigenValues <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/joint_gal_mainland_biallelic_snp_pca_results.eigenval", delim = " ", col_names = F)

eigenVectors <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/joint_gal_mainland_biallelic_snp_pca_results.eigenvec", delim = "\t", col_names = T) 
# pop_map.txt: at least two columns: IID, population
pop_map <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_pop_data.csv") 

#env_data <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/sample_population_info.csv")

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

eigen_env <- left_join(pop_map, eigenVectors, by = "IID")

my_palette <- c(
  Outlier = "#E41A1C",
  Mainland = "#377EB8",
  Northern = "#4DAF4A",
  Western = "#984EA3",
  Central_SE = "#FF7F00",
  `Far North` = "#8c8c8c"
)

#Plot PCA
ggplot(data = eigen_env,mapping = aes(x = PC1, y = PC2, color = Bioregion))+#, shape = Cluster)) +
  geom_point(size = 3, show.legend = T,alpha=0.8 ) +
  scale_color_manual(values = my_palette) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Gal & Mainland (SNPs 29271149)",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()


#######PCA No OUTLIER##########
eigenValues <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/joint_biallelic_snp_pca_No_Outliers.eigenval", delim = " ", col_names = F)

eigenVectors <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/joint_biallelic_snp_pca_No_Outliers.eigenvec", delim = "\t", col_names = T) 

pop_map <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_pop_data.csv") |> 
  filter(!IID %in% c("GAL-225","GAL-057","GAL-127")) #remove duplicate sample



## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

eigen_env <- left_join(pop_map, eigenVectors, by = "IID") 

ggplot(data = eigen_env,mapping = aes(x = PC1, y = PC2, color = Bioregion))+#, shape = Cluster)) +
  geom_point(size = 3, show.legend = T,alpha=0.8 ) +
  scale_color_manual(values = my_palette) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Gal & Mainland (No Outliers - SNPs 28822616)",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()

##HET

# Load heterozygosity
het <- read_table("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/het_joint_qc_No_Outliers.het", col_types = cols())

het <- het |>
  
  mutate(
    obs_het = 1 - (`O(HOM)` / OBS_CT),   # observed heterozygosity
    exp_het = 1 - (`E(HOM)` / OBS_CT)    # expected heterozygosity
  ) 

pop_map <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_pop_data.csv") |> 
  filter(!IID %in% c("GAL-225","GAL-057","GAL-127")) #remove duplicate sample



het_pop <- het |>
  left_join(pop_map, by = "IID")

### 3. Boxplot of observed heterozygosity by population
ggplot(het_pop, aes(x = Cluster, y = obs_het, fill = Cluster)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Population",
    y = "Observed heterozygosity",
    title = "Per-sample observed heterozygosity by cluster (No Outliers)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### 5. Optional: plot inbreeding coefficient F by population
ggplot(het_pop, aes(x = Cluster, y = F, fill = Cluster)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Population",
    y = "Inbreeding coefficient (F)",
    title = "Per-sample F by cluster (No Outliers)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

########################
#######################
####RIGHT CLUSTER ONLY PCA####

eigenValues <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/joint_RIGHT_snp_pca_No_Outliers.eigenval", delim = " ", col_names = F)
eigenVectors <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/joint_Genotype_final/joint_RIGHT_snp_pca_No_Outliers.eigenvec", delim = "\t", col_names = T)
pop_map <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_pop_data.csv") |> 
  filter(Cluster %in% c("RIGHT","Mainland")) #filter RIGHT CLUSTER ONLY
## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)
eigen_env <- left_join(pop_map, eigenVectors, by = "IID")
#Plot PCA
ggplot(data = eigen_env,mapping = aes(x = PC1, y = PC2, color = Bioregion))+#, shape = Cluster)) +
  geom_point(size = 3, show.legend = T ,alpha=0.8 ) +
  scale_color_manual(values = my_palette) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Gal & Mainland - RIGHT Cluster (38 samples)",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()


















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



