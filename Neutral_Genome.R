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

eigenValues <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/Neutral_Gal_Only/full_neutral_genome_snp_plink_biallelic_snp_pca.eigenval", delim = " ", col_names = F) 

eigenVectors <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/Neutral_Gal_Only/full_neutral_genome_snp_plink_biallelic_snp_pca.eigenvec", delim = "\t", col_names = T) 

env_data <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_pop_data.csv")

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

eigen_env <- left_join(eigenVectors, env_data, by = "IID") |> 
  filter(!IID %in% c("GAL-057","GAL-127")) #|> 
  # mutate(clust=case_when(PC2 > -0.05 ~ "RIGHT",
  #                        PC2 < -0.10 ~ "LEFT",
  #                        PC2 < -0.05 & PC2 > -0.10 ~ "CENTER"))

#write_csv(eigen_env,"/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/combined_gal_mainland_no_57_127_pca_env_data.csv")

library(RColorBrewer)
populations <- c('Central_SE', 'Far North', 'Mainland',  'Outlier', 'Northern', 'Western')
my_palette <- setNames(brewer.pal(6, "Dark2"), populations)
my_palette <- c(
  Outlier = "#E41A1C",
  Mainland = "#377EB8",
  Northern = "#4DAF4A",
  Western = "#984EA3",
  Central_SE = "#FF7F00",
  `Far North` = "#8c8c8c"
)

my_palette2 <- c(
  RIGHT = "#53051D",
  Center = "#F12761",
  LEFT = "#F2CE1B"
)

#Plot PCA
ggplot(data = eigen_env,mapping = aes(x = -PC1, y = PC2, color = new_Cluster)) +
  geom_point(size = 3, show.legend = T,alpha=0.8 ) +
  scale_color_manual(values = my_palette2) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Neutral Gal- No Outliers (57 & 127) - SNPs 23634852",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()


########## LD Pruned###############
eigenValues <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/Neutral_Gal_Only/final_pca_resultsLD_20_5_0p1.eigenval", delim = " ", col_names = F)

eigenVectors <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/Neutral_Gal_Only/final_pca_resultsLD_20_5_0p1.eigenvec", delim = "\t", col_names = T)

env_data <- read_csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/Neutral_Gal_Only/sample_population_info.csv")

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

eigen_env <- left_join(eigenVectors, env_data, by = "IID") |> 
  filter(!IID %in% c("GAL-057","GAL-127"))

#Plot PCA
ggplot(data = eigen_env,mapping = aes(x = PC1, y = PC2, color = Cluster, shape = Bioregion)) +
  geom_point(size = 3, show.legend = T ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Neutral Gal - LD Pruned - SNPs 6411176",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()




########## PCA per Global cluster (left, right, center)

left_cluster <- eigen_env |> 
  filter(Cluster=="LEFT")

ggplot(data = left_cluster,mapping = aes(x = PC1, y = PC2, color = Bioregion, shape = Bioregion)) +
  geom_point(size = 3, show.legend = T ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Neutral Gal - LEFT Cluster",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()

#Right cluster
right_cluster <- eigen_env |> 
  filter(Cluster=="RIGHT")

ggplot(data = right_cluster,mapping = aes(x = PC1, y = PC2, color = Bioregion, shape = Bioregion)) +
  geom_point(size = 3, show.legend = T ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") +
  #scale_color_viridis_c() +
  labs(title = "PCA Neutral Gal - RIGHT Cluster",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()
#Center cluster
center_cluster <- eigen_env |> 
  filter(Cluster=="Center")

ggplot(data = center_cluster,mapping = aes(x = PC1, y = PC2, color = Bioregion, shape = Bioregion)) +
  geom_point(size = 3, show.legend = T ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  #geom_text(hjust=0, vjust=0,aes(label=IID)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(legend.position="right") + 
  #scale_color_viridis_c() +
  labs(title = "PCA Neutral Gal - CENTER Cluster",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"))+
  theme_minimal()
