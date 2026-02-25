library(sf)

# read las
las <- readLAS("PC/Wenns_small.las")

# get projection
wkt <- paste(readLines("PC/WennsUTM.prj"), collapse = "")
crs <- st_crs(wkt)

# define projection of las
st_crs(las) <- crs

aoi <- read_sf("AOI/AOI.shp")

las_aoi <- clip_roi(las, aoi)

writeLAS(las_aoi, "PC/Wenns_aoi.las")