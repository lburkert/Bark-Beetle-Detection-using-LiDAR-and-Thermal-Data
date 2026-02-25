################################################################
#                                                              #
# SCRIPT TO Perform Individual Tree Detection and Segmentation #
#                                                              #
################################################################

# Clear environment
rm(list=ls())

# Libraries 
library(this.path)
library(lidR)
library(sf)
library(ggplot2)


# Set path
script_path=this.path()
setwd(dirname(script_path))

####
# Insert Crop_aoi.r
####

# read las file
las_aoi <- readLAS("PC/Wenns_aoi.las")     

#### Statistical outlier removal ####
# apply SOF to las
las_nonoise = classify_noise(las_aoi, sor(15,7))

# get only noise points and save them
las_noise <- filter_poi(las_nonoise, Classification == 18)
writeLAS(las_noise, "noise_points.las")

# ---- get cleaned point cloud ----
las_clean <- filter_poi(las_nonoise, Classification != 18)



####  Cloth Simulation Filter ####
# classify ground
las_ground <- classify_ground(las_clean, algorithm = csf(cloth_resolution = 0.1))
las_ground <- classify_ground(las_clean, algorithm = pmf())

p1 <- c(las_ground@header@PHB$`Min X`, las_ground@header@PHB$`Min Y`)
p2 <- c(las_ground@header@PHB$`Max X`, las_ground@header@PHB$`Max Y`)
las_tr <- clip_transect(las_ground, p1, p2, width = 0.5, xz = TRUE)

x_range <- range(payload(las_tr)$X)
z_range <- range(payload(las_tr)$Z)

# Plot
p <- ggplot(payload(las_tr), aes(X, Z, color = factor(Classification))) +
  geom_point(size = 0.15) +
  coord_equal(xlim = x_range, ylim = z_range, expand = FALSE) +  # tight cropping
  scale_color_manual(
    values = c(
      "0" = "#1b7837",   # Non-ground
      "2" = "#8c510a"    # Ground
    ),
    labels = c(
      "0" = "Non-ground",
      "2" = "Ground"
    )
  ) +
  labs(
    x = "Distance (m)",
    y = "Elevation (m)",
    color = "Classification"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4, shape = 16) 
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.key = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "pt")
  )


ggsave(
  filename = "plots/las_csf_01.png",
  plot = p,
  width = 30,
  height = 9,
  dpi = 700,       
  units = "cm",
  bg = "white"
)
  

gnd <- filter_ground(las_ground)
plot(gnd, size = 3, bg = "white") 
plot(las_ground, color = "Classification", size = 3, bg = "white")

# Digital terrain model
dtm_tin <- rasterize_terrain(las_ground, res = 0.2, algorithm = tin())
plot_dtm3d(dtm_tin, bg = "white") 

dtm_idw <- rasterize_terrain(las_ground, res = 0.2, knnidw())
plot_dtm3d(dtm_idw, bg = "white")

# height normalization
nlas <- las_ground - dtm_idw
plot(nlas, size = 4, bg = "white")
writeLAS(nlas, "nomalized_aoi.las")

# check ground points
hist(filter_ground(nlas)$Z, breaks = seq(-1.5, 1.5, 0.01), main = "", xlab = "Elevation")

#### Calculate Canopy Height Model ####
chm <- rasterize_canopy(nlas, 0.2)
col <- height.colors(25)
plot(chm, col = col)

# Individual Tree Detection
ttops_5 <- locate_trees(las_ground, lmf(ws = 5))
ttops_2 <- locate_trees(las_ground, lmf(ws = 2))
ttops_1 <- locate_trees(las_ground, lmf(ws = 1))
ttops_25 <- locate_trees(las_ground, lmf(ws = 2.5))

plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops_5), add = TRUE, pch = 3)





# Point-to-raster 2 resolutions
chm_p2r_02 <- rasterize_canopy(nlas, 0.2, p2r(subcircle = 0.2), pkg = "terra")
chm_p2r_1 <- rasterize_canopy(nlas, 1, p2r(subcircle = 0.2), pkg = "terra")

# Pitfree with and without subcircle tweak
chm_pitfree_05_1 <- rasterize_canopy(nlas, 0.5, pitfree(), pkg = "terra")
chm_pitfree_05_2 <- rasterize_canopy(nlas, 0.5, pitfree(subcircle = 0.2), pkg = "terra")
chm_pitfree_05_3 <- rasterize_canopy(nlas, 0.2, pitfree(), pkg = "terra")

# Post-processing median filter
kernel <- matrix(1,3,3)
chm_p2r_02_smoothed <- terra::focal(chm_p2r_02, w = kernel, fun = median, na.rm = TRUE)
chm_p2r_1_smoothed <- terra::focal(chm_p2r_1, w = kernel, fun = median, na.rm = TRUE)

ttops_chm_p2r_02 <- locate_trees(chm_p2r_02, lmf(2))
ttops_chm_p2r_1 <- locate_trees(chm_p2r_1, lmf(5))
ttops_chm_pitfree_05_1 <- locate_trees(chm_pitfree_05_1, lmf(5))
ttops_chm_pitfree_05_2 <- locate_trees(chm_pitfree_05_2, lmf(5))
ttops_chm_pitfree_05_3 <- locate_trees(chm_pitfree_05_3, lmf(5))
ttops_chm_p2r_02_smoothed <- locate_trees(chm_p2r_02_smoothed, lmf(2.5))
ttops_chm_p2r_1_smoothed <- locate_trees(chm_p2r_1_smoothed, lmf(5))

par(mfrow=c(3,2))
col <- height.colors(50)
plot(chm_p2r_02, main = "CHM P2R 0.5", col = col); plot(sf::st_geometry(ttops_chm_p2r_02), add = T, pch =3)
plot(chm_p2r_1, main = "CHM P2R 1", col = col); plot(sf::st_geometry(ttops_chm_p2r_1), add = T, pch = 3)
plot(chm_p2r_02_smoothed, main = "CHM P2R 0.5 smoothed", col = col); plot(sf::st_geometry(ttops_chm_p2r_02_smoothed), add = T, pch =3)
plot(chm_p2r_1_smoothed, main = "CHM P2R 1 smoothed", col = col); plot(sf::st_geometry(ttops_chm_p2r_1_smoothed), add = T, pch =3)
plot(chm_pitfree_05_1, main = "CHM PITFREE 1", col = col); plot(sf::st_geometry(ttops_chm_pitfree_05_1), add = T, pch =3)
plot(chm_pitfree_05_2, main = "CHM PITFREE 2", col = col); plot(sf::st_geometry(ttops_chm_pitfree_05_2), add = T, pch =3)


x <- plot(nlas, col = rgb, bg = "white", size = 4)
add_treetops3d(x, ttops_chm_p2r_02_smoothed)

# XYZ extrahieren
coords <- st_coordinates(ttops_chm_p2r_02_smoothed)

df <- data.frame(
  X = coords[,1],
  Y = coords[,2],
  Z = coords[,3]
)

las <- LAS(df)

projection(las) <- st_crs(ttops_chm_p2r_02_smoothed)$wkt

writeLAS(las, "ttops_2.5.las")

las_check(las_aoi)


# ---- Individual Tree Segmentation ----

algo <- dalponte2016(chm_p2r_02_smoothed, ttops_chm_p2r_02_smoothed)
las <- segment_trees(nlas, algo) # segment point cloud
plot(las, bg = "white", size = 4, color = "treeID") # visualize trees
writeLAS(las, "tree_segmented.las")
