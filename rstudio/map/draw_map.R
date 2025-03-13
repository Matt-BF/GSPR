library("tidyverse")
library("dbscan")
library("sf")
library("geosphere")
library("cluster")
library("geosphere")
library("rnaturalearth")
library("rnaturalearthdata")


coordinates <- read.csv("df_map_fixed.csv")

unique_coordinates <- coordinates %>%
  select(Latitude, Longitude) %>%
  unique()

# Convert to spatial points dataframe
coords_sf <- st_as_sf(unique_coordinates, coords = c("Longitude", "Latitude"), crs = 4326)

# Compute distance matrix in meters (using geosphere::distm)
dist_matrix <- distm(st_coordinates(coords_sf))

# Apply DBSCAN clustering
eps <- 5e+6 # Set epsilon (distance cutoff) in meters
db <- dbscan(dist_matrix, eps = eps, minPts = 1)

# Add cluster labels to original data
unique_coordinates["cluster"] <- db$cluster

# Join the dataframes to put the cluster information in the one containing votus
coordinates <- left_join(coordinates, unique_coordinates, by = c("Latitude", "Longitude"))

# Compute the number of distinct votus per cluster
n_ptus <- coordinates %>%
  group_by(cluster) %>%
  summarise(n = n_distinct(PTU))

# Function to find medoid
find_medoid <- function(cluster_points) {
  if (nrow(cluster_points) == 1) {
    return(cluster_points)
  }
  dist_matrix_cluster <- distm(st_coordinates(st_as_sf(cluster_points, coords = c("Longitude", "Latitude"), crs = 4326)))
  medoid_index <- pam(dist_matrix_cluster, k = 1)$id.med
  return(cluster_points[medoid_index, , drop = FALSE])
}

# Apply the function to each cluster
medoids <- coordinates %>%
  select(Latitude, Longitude, cluster) %>%
  unique() %>%
  group_by(cluster) %>%
  group_modify(~ find_medoid(.x)) %>%
  ungroup()

# Join dataframes to have the number of distinct ptus together with the medoids
medoids <- left_join(medoids, n_ptus, by = "cluster")

aeqd_proj <- st_crs(paste0("+proj=loxim"))

# Transform medoid points to the new projection
medoids_sf <- st_as_sf(medoids, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(aeqd_proj) %>%
  mutate(cluster = as_factor(cluster))

# Load world map data
world <- ne_countries(scale = "small", returnclass = "sf") %>%
  st_transform(aeqd_proj)


# Plot the map with medoid points
p1 <- ggplot() +
  geom_sf(data = world, fill = "#EFE1D9", color = NA) +
  geom_sf(data = medoids_sf, aes(size = n), color = "#277F8E", alpha = 0.5, stroke = 0) +
  #scale_size_continuous(breaks = seq(1000, 6000, by=2000))+
  theme_minimal() +
  theme(
    legend.position = "bottom",
  )
print(p1)
ggsave("ptu_map.pdf", height = 250, width = 420, units = "px", dpi = 72)
