# Load packages
library(sf)
library(ggplot2)
library(dplyr)
library(ggrepel)

#setwd("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER")

# Read the CSV file
data <- read.csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_pop_data.csv",header=T) |> 
  filter(!IID %in% c("GAL-057","GAL-127","GAL-225")) #remove outliers

data_clean <- na.omit(data)
data_clean <- data_clean[data_clean$isla!="Salango",]


# Add a 'Count' variable and maintain all other columns
data_with_count <- data_clean %>%
  group_by(clust,GPS_lat,GPS_lon) %>%
  mutate(Count = n()) %>%
  ungroup() %>%
  group_by(clust,GPS_lat,GPS_lon, .add = TRUE) %>%
  slice(1) %>%
  ungroup()


# Convert the summarized data to an sf object
points_sf <- st_as_sf(data_with_count, coords = c("GPS_lon", "GPS_lat"), crs = 4326)

#Obtained Galaapgos shape file from: https://www.protectedplanet.net/187
# Read the shapefile (replace 'path_to_shapefile' with the actual path to your .shp file)
shp_file <- st_read("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/WDPA_WDOECM_Mar2024_Public_187_shp_0/WDPA_WDOECM_Mar2024_Public_187_shp-polygons.shp")


# Convert the summarized data to a regular data frame for plotting with ggrepel
points_df <- as.data.frame(data_with_count)

my_palette <- c(
  RIGHT = "#53051D",
  Center = "#F12761",
  LEFT = "#F2CE1B"
)

# Plot the shapefile as the base layer and the points on top
map_plot <- ggplot() +
  geom_sf(data = shp_file) +
  geom_point(data = points_df, aes(x = GPS_lon, y = GPS_lat, color = new_Cluster),position = position_jitter(width = 0.05, height = 0, seed = 123), alpha = 0.85, show.legend = T,size=2) +
  scale_color_manual(values = my_palette) +
  #ggrepel::geom_text_repel(data = points_df, aes(x = GPS_lon, y = GPS_lat), size = 2) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "grey", size = 0.25, linetype = "dotted"),
        panel.grid.minor = element_line(color = "grey", size = 0.25, linetype = "dotted"),
        legend.position = "right")


map_plot


##### Plot the shapefile as the base layer and the points on top #########
ggplot() +
  geom_sf(data = shp_file) +
  geom_sf(data = points_sf, aes(size = Count, color = clust), alpha = 0.6) +
  geom_sf_text(data = points_sf, aes(label = Count), hjust = 1, vjust = 1, color = "black") +
  scale_size_continuous(range = c(3, 12)) +
  labs(x = "Longitude", y = "Latitude", size = "Count") +
  theme_minimal() +
  theme(legend.position = "none") # Hide the legend if not needed



########## Ecuador map bb(-86.13472477682836,-9.145419142341918,-68.5111865619435,3.3653202060978) #####
shp_ec <- st_read("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/ec_shp/ec.shp")

ggplot() +
  geom_sf(data = shp_ec)

ec_crop=st_crop(shp_ec, xmin = -86.13472477682836, ymin = -9.145419142341918, xmax = -68.5111865619435, ymax = 3.3653202060978)

ggplot() +
  geom_sf(data = ec_crop)

##South America /Users/jaimeortiz/Downloads/stanford-vc965bq8111-shapefile/vc965bq8111.shp
shp_ec <- st_read("/Users/jaimeortiz/Downloads/stanford-vc965bq8111-shapefile/vc965bq8111.shp")

ggplot() +
  geom_sf(data = shp_ec)
#Make geometries valid
shp_ec_valid <- st_make_valid(shp_ec)

# Now crop
ec_crop <- st_crop(shp_ec_valid, xmin = -82.4604, ymin = -59.48714, xmax = -35.23419, ymax = 12.62908)

ggplot() +
  geom_sf(data = ec_crop)

# Read the CSV file
samples <- read.csv("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/new_combined_gal_mainland_no_57_127_pca_env_data.csv",header=T)

df_clean <- na.omit(samples)

# Add a 'Count' variable and maintain all other columns
data_with_count <- df_clean %>%
  group_by(clust,GPS_lat,GPS_lon) %>%
  mutate(Count = n()) %>%
  ungroup() %>%
  group_by(clust,GPS_lat,GPS_lon, .add = TRUE) %>%
  slice(1) %>%
  ungroup()


# Convert the summarized data to an sf object
points_sf <- st_as_sf(data_with_count, coords = c("GPS_lon", "GPS_lat"), crs = 4326)

# Convert the summarized data to a regular data frame for plotting with ggrepel
points_df <- as.data.frame(data_with_count)

# Plot the shapefile as the base layer and the points on top
salango_plot <- ggplot() +
  geom_sf(data = ec_crop) +
  geom_point(data = points_df[points_df$isla=="Salango",], aes(x = GPS_lon, y = GPS_lat), size=5,color = "lightblue", alpha = 0.8, show.legend = FALSE) +
  #ggrepel::geom_text_repel(data = points_df[points_df$isla=="Salango",], aes(x = GPS_lon, y = GPS_lat, label = Count), size = 4) +
  #labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  #theme(panel.grid.major = element_line(color = "grey", size = 0.25, linetype = "dotted"),
  #panel.grid.minor = element_line(color = "grey", size = 0.25, linetype = "dotted")) +
  xlim(c(-85, NA)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


salango_plot

#### combine plots ####

#install.packages("cowplot")
# Load the packages
library(cowplot)

# Add a border to the inset plot by modifying its theme
salango_plot <- salango_plot + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

salango_plot <- salango_plot +
  theme(
    #axis.text.x = element_blank(),  # Remove X axis text
    #axis.text.y = element_blank(),  # Remove Y axis text
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(), # Remove X axis title
    axis.title.y = element_blank() # Remove Y axis title
    #axis.ticks = element_blank()    # Remove axis ticks
  )

# Create the combined plot with the inset in the upper right corner
combined_plot <- ggdraw(map_plot) +
  draw_plot(salango_plot, x = 0.4, y = 0.45, width = 0.9, height = 0.5) # Adjust these values as needed

# Display the combined plot
print(combined_plot)

