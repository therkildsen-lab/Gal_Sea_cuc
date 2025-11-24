#load libraries
library(tidyverse) |> suppressPackageStartupMessages()
library(readr)
library(tidyr)
library(dplyr)
#devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
library(data.table)
library(stringr)
library(viridis)



#Filtered VCF file No samples 57 and 127
#read eigenvectors and eigenvalues (Produced using PLINK2)

eigenValues <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/DP100_mainland_biallelic_snp_pca_results.eigenval", delim = " ", col_names = F)

eigenVectors <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/bcftools_merged_PLINK/DP100_mainland_biallelic_snp_pca_results.eigenvec", delim = "\t", col_names = T) |> 
  rename(IID = `#IID`)

env_data <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_no_57_127_pca_env_data.csv")

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

eigen_env <- left_join(eigenVectors, env_data, by = "IID") |> 
  filter(!IID %in% c("GAL-225"),clust!="LEFT") 
# |> 
#   mutate(clust=case_when(PC2 > -0.05 ~ "RIGHT",
#                          PC2 < -0.10 ~ "LEFT",
#                          PC2 < -0.05 & PC2 > -0.10 ~ "CENTER"))

#write_csv(eigen_env,"/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/combined_gal_mainland_no_57_127_pca_env_data.csv")


#Plot PCA
ggplot(data = eigen_env,mapping = aes(x = PC1, y = PC2, color = Bioregion, shape = clust)) +
  geom_point(size = 3, show.legend = T ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Right Cluster",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()
