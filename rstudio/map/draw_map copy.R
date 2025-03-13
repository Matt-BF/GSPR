library("tidyverse")
library("dbscan")
library("sf")
library("geosphere")
library("cluster")
library("geosphere")
library("rnaturalearth")
library("rnaturalearthdata")

coordinates <- read_tsv("uhgv_metadata.tsv", na = c("", "NULL")) %>%
  mutate(
    latitude = if_else(latitude == 0 & longitude == 0, NA, latitude),
    longitude = if_else(latitude == 0 & longitude == 0, NA, longitude),
    latitude = if_else(latitude == 33.86 & longitude == 151.2111, NA, latitude),
    longitude = if_else(latitude == 33.86 & longitude == 151.2111, NA, longitude),
    latitude = if_else(is.na(latitude), country_latitude, latitude),
    longitude = if_else(is.na(longitude), country_longitude, longitude)
  ) %>%
  select(uhgv_votu, latitude, longitude) %>%
  drop_na()

unique_coordinates <- coordinates %>%
  select(latitude, longitude) %>%
  unique()

# Convert to spatial points dataframe
coords_sf <- st_as_sf(unique_coordinates, coords = c("longitude", "latitude"), crs = 4326)

# Compute distance matrix in meters (using geosphere::distm)
dist_matrix <- distm(st_coordinates(coords_sf))

# Apply DBSCAN clustering
eps <- 4200000 # Set epsilon (distance cutoff) in meters
db <- dbscan(dist_matrix, eps = eps, minPts = 1)

# Add cluster labels to original data
unique_coordinates["cluster"] <- db$cluster

# Join the dataframes to put the cluster information in the one containing votus
coordinates <- left_join(coordinates, unique_coordinates, by = c("latitude", "longitude"))

# Compute the number of distinct votus per cluster
n_votus <- coordinates %>%
  group_by(cluster) %>%
  summarise(n = n_distinct(uhgv_votu))

# Function to find medoid
find_medoid <- function(cluster_points) {
  if (nrow(cluster_points) == 1) {
    return(cluster_points)
  }
  dist_matrix_cluster <- distm(st_coordinates(st_as_sf(cluster_points, coords = c("longitude", "latitude"), crs = 4326)))
  medoid_index <- pam(dist_matrix_cluster, k = 1)$id.med
  return(cluster_points[medoid_index, , drop = FALSE])
}

# Apply the function to each cluster
medoids <- coordinates %>%
  select(latitude, longitude, cluster) %>%
  unique() %>%
  group_by(cluster) %>%
  group_modify(~ find_medoid(.x)) %>%
  ungroup()

# Join dataframes to have the number of distinct votus together with the medoids
medoids <- left_join(medoids, n_votus, by = "cluster")

# Define the azimuthal equidistant projection centered on Brazil
aeqd_proj <- st_crs(paste0("+proj=aeqd +lat_0=", -7.0, " +lon_0=", -20.0))

# Transform medoid points to the new projection
medoids_sf <- st_as_sf(medoids, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(aeqd_proj) %>%
  mutate(cluster = as_factor(cluster))

# Load world map data
world <- ne_countries(scale = "small", returnclass = "sf") %>%
  st_transform(aeqd_proj)

# Define the circle for the globe stroke with explicit CRS
globe_circle <- st_sfc(
  st_buffer(
    st_point(c(0, 0)),
    dist = max(abs(st_bbox(world)))
  ) * 1.07,
  crs = aeqd_proj
)

# Plot the map with medoid points
ggplot() +
  geom_sf(data = globe_circle, fill = "#E8F6FC", color = NA) +
  geom_sf(data = world, fill = "#E5D4AE", color = NA) +
  geom_sf(data = medoids_sf, aes(size = n), color = "#E25A63", alpha = 0.5, stroke = 0) +
  scale_size_continuous(breaks = c(1000, 10000, 20000, 40000)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Geographical distribution of virus genomes",
    subtitle = "Azimuthal Equidistant Projection",
    size = "No. distinct vOTUs"
  )
ggsave("map1.pdf", height = 250, width = 420, units = "px", dpi = 72)
