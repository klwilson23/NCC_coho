library(dbplyr)

buffer <- 25000
library(bcdata)
library(dplyr)
library(bcmaps)
library(sf)
library(sp)
library(ggplot2)
library(ggspatial) # this is for adding the arrows
library(rgdal) #use this for data conversion
library(ggrepel) # to offset labels using geom_sf_label_repel  --> Not done here
library(riverdist) # to snap points to River --> Not done here
library(bcmapsdata)
library(viridis)
library(ggnewscale)
library(tidyr)

NCC_coho <- read.csv("Data/NCC coho.csv",header=TRUE)

NCC_coho$lon <- NCC_coho$Longitude
NCC_coho$lat <- NCC_coho$Latitude

data2 <- NCC_coho[,c("lon","lat")] #%>% mutate(UTM_W=as.numeric(UTM_W),UTM_N=as.numeric(UTM_N)) #
sputm <- SpatialPoints(data2, proj4string=CRS("+proj=longlat"))
spgeo <- spTransform(sputm, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
data_points <- as_tibble(NCC_coho) %>% mutate(lat=spgeo@coords[,2],lon=spgeo@coords[,1]) %>%
  st_as_sf(coords=c("lon","lat"),crs=3005)

## Subset the points data to be a single point for each location (for labelling Locations)
data_point_labels <- data_points %>% group_by(population) %>% slice(1)

ncc_extent <- raster::extent(data_points)
ncc_extent@xmin <- ncc_extent@xmin - 3*buffer
ncc_extent@ymin <- ncc_extent@ymin - buffer 
ncc_extent@xmax <- ncc_extent@xmax + buffer 
ncc_extent@ymax <- ncc_extent@ymax + buffer 
plot_area_ncc <- ncc_extent %>%
  st_bbox() %>%                 # Turn into a square
  st_as_sfc(crs=st_crs(data_points))

st_crs(plot_area_ncc) <- st_crs(data_points)

watersheds <- bcdc_query_geodata('51f20b1a-ab75-42de-809d-bf415a0f9c62') %>%
  filter(INTERSECTS(plot_area_ncc)) %>%      # not sure about this line
  collect() %>%                           #Extracts the data
  st_intersection(plot_area_ncc)             #Where it intersects with plot line


SR_dat <- readRDS(file="Data/time_varying_ocean_entry.rds")
log_alpha <- read.csv("Data/Coho alpha trends.csv")
SR_datsub <- SR_dat[,c(1:30,37)]

SR_regional <- SR_datsub %>%
  group_by(year,region) %>%
  summarise(geometry=mean(geometry))

df_full <-merge(SR_regional,log_alpha,by.x=c("year","region"),by.y=c("year","region"))
df_full$x <- st_coordinates(df_full)[,1]
df_full$y <- st_coordinates(df_full)[,2]
df_full$crs <- st_crs(df_full)[[1]]
write.csv(df_full,"Data/time-varying ocean entry.csv",row.names=FALSE)
