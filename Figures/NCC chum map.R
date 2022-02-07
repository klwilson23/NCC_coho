# See: https://github.com/bcgov/bcmaps
# see: https://cran.r-project.org/web/packages/bcmaps/bcmaps.pdf
# availble layers in bcmaps:
# https://gist.github.com/ateucher/86674e12ffba66ecce87cceb7bffbf41
# https://github.com/poissonconsulting/fwabc#<https://github.com/poissonconsulting/fwabc>
# https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html

buffer <- 15000

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
# Set plot box.  Set Port Hardy as centre and box is 20km to each side
plot_area1 <- bc_cities() %>%   #Extract datafram of all bc cities location
  filter(NAME == "Bella Bella") %>%
  #filter(NAME == "Port Hardy") %>%
  st_buffer(dist = 1000000)%>%    # 20000 is a decent zoom in one that covers the Keogh
  st_bbox() %>%                 # Turn into a square
  st_as_sfc()                   # converts to sfc object
NCC_chum <- read.csv("Data/NCC chum.csv",header=TRUE)

NCC_chum$lon <- NCC_chum$long
data2 <- NCC_chum[,c("lon","lat")] #%>% mutate(UTM_W=as.numeric(UTM_W),UTM_N=as.numeric(UTM_N)) #
sputm <- SpatialPoints(data2, proj4string=CRS("+proj=longlat"))
spgeo <- spTransform(sputm, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
data_points <- as_tibble(NCC_chum) %>% mutate(lat=spgeo@coords[,2],lon=spgeo@coords[,1]) %>%
  st_as_sf(coords=c("lon","lat"),crs=3005)

## Subset the points data to be a single point for each location (for labelling Locations)
data_point_labels <- data_points %>% group_by(population) %>% slice(1)

ncc_extent <- raster::extent(data_points)
ncc_extent@xmin <- ncc_extent@xmin - buffer
ncc_extent@ymin <- ncc_extent@ymin - buffer 
ncc_extent@xmax <- ncc_extent@xmax + buffer 
ncc_extent@ymax <- ncc_extent@ymax + buffer 
plot_area_ncc <- ncc_extent %>%
  st_bbox() %>%                 # Turn into a square
  st_as_sfc(crs=st_crs(data_points))
st_crs(plot_area_ncc) <- st_crs(data_points)
#bcdc_get_record("https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-rivers")
rivers_in_plot_area <- bcdc_query_geodata('f7dac054-efbf-402f-ab62-6fc4b32a619e') %>%
  #filter(STREAM_ORDER %in% c(3,4,5)) %>%  #Defines as only streams order 3,4,5 (too many including 1 &2)
  filter(INTERSECTS(plot_area_ncc)) %>%      # not sure about this line
  collect() %>%                           #Extracts the data
  st_intersection(plot_area_ncc)             #Where it intersects with plot line

rivernames <- NCC_chum$population
rivers_ncc <- rivers_in_plot_area %>% filter(GNIS_NAME_1 %in% rivernames)

ncc_wshed <- rivers_ncc %>%
  filter(WATERSHED_GROUP_ID%in%unique(rivers_ncc$WATERSHED_GROUP_ID))  %>% 
  filter(LEFT_RIGHT_TRIBUTARY!="NONE")

## --- Collect Mapping layers
# below grabs the coastline data within plot box
coast_line <- bc_bound_hres() %>% st_intersection(plot_area_ncc)

## coastlines taken from freshwater atlas
## This looks the same as above but takes a lot longer to run
#coastline_k <- bcdc_query_geodata("freshwater-atlas-coastlines") %>% collect() %>%  st_intersection(square_round_ph)

# Locations of all cities within plot box
city_df<- bc_cities() %>% st_intersection(plot_area_ncc)

#Below is to get a dataset of lakes in plot box
lakes_in_plot_area <- bcdc_query_geodata("freshwater-atlas-lakes") %>%
  filter(INTERSECTS(plot_area_ncc)) %>%
  collect() %>%
  st_intersection(plot_area_ncc)

# Ocean colouring - this doesn't work well  because the resolution isn't the same
ocean_colour <- bc_neighbours() %>% 
  filter(type=="Ocean") %>% 
  st_intersection(plot_area_ncc)
###########

# load forestry data
library(raster)
str_name<-'~/Google Drive/SFU postdoc/Keogh river/BC_disturbance/logging_ageclass2012/logging_year.tif' 
log_year <- raster(str_name)
# project forestry data into map projections
proj4string(log_year) <- CRS("+init=EPSG:3005")
# find the map extents from the raster file, and "crop" or subset the raster, otherwise its a huge file
ext <- extent(plot_area_ncc[[1]][1][[1]][1:4,])
log_year_crop <- crop(log_year,ext)
log_year_crop@data@attributes[[1]]$Rowid <- log_year_crop@data@attributes[[1]]$ID
logging_df <- as.data.frame(log_year_crop,xy=TRUE)
colnames(logging_df) <- c("x", "y","Year","count")
logging_df$Year <- as.numeric(logging_df$Year)
logging_df <- logging_df[complete.cases(logging_df),]
# maniuplate logging data to bin last logged year by 10 year bins
logging_df$Year <- round(logging_df$Year/10)*10
# anything earlier than 1950 was last logged "Pre-1950"
logging_df$Year <- ifelse(logging_df$Year<=1950,"Pre-1950",logging_df$Year)
logging_df$Year <- factor(logging_df$Year,levels=c("Pre-1950",sort(unique(logging_df$Year)[!grepl("Pre",unique(logging_df$Year))])))

## The fun part! we get to make a map
## Plot map -- Order of plotting MATTERS
ncc_map <- ggplot() +
  geom_sf(data = coast_line, fill = "grey90") +                 #Plot coastline
  geom_sf(data = plot_area_ncc, alpha = 0,colour='black') +        #Plot area box
  geom_tile(data=logging_df, aes(x=x, y=y, fill=Year), alpha=0.8) +
  #scale_fill_viridis(name = "Last logging year",option="viridis",direction=-1) +
  #scale_fill_gradient2(name="Last logged year",midpoint = 1920, low="darkgreen",mid="olivedrab3",high="goldenrod1",space="Lab") +
  scale_fill_manual(name="Last logged year",values=c("#1b7837","#7fbf7b","#01665e","#5ab4ac","#762a83","#d8b365","#8c510a")) +
  #scale_fill_brewer(name="Last logged year",palette = "Paired",direction=-1) +
  geom_sf(data = rivers_in_plot_area, colour = "lightblue3") +  #Plot Rivers
  geom_sf(data = rivers_ncc, colour = "lightblue3") +
  #geom_sf(data = lakes_in_plot_area, fill = "lightblue1") +     #Plot Lakes
  ggsflabel::geom_sf_label_repel(data=data_point_labels,aes(label=population),size=1.5, force = 1, nudge_x = -2, seed = 10)+ #Add labels
  coord_sf(expand = FALSE) +                                    #Expands box to axes
  xlab('Longitude') + ylab('Latitude') +                        #Axis titles
  annotation_scale(location = "br", width_hint = 0.5) +         #Rose Compass
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.5, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(1,"cm"), width = unit(1, "cm"))+
  theme(panel.background = element_rect('lightblue1')           #Make empty space blue to colour ocean
        , panel.grid.major = element_line('lightblue1'),
        legend.position="top",legend.text = element_text(size=6),legend.key.width=unit(0.35, "in"))
bec <- bec() %>% st_intersection(plot_area_ncc) %>% st_intersection(bc_bound_hres())
#bec <- bec[! (bec$ZONE %in% "CDF"),]
pfma <- st_read("Data/PFMA/DFO_PFMA_SUBAREAS_SP/DFO_SBAREA_polygon.shp")
pfma <- pfma[pfma$MGMT_AREA>=6 & pfma$MGMT_AREA<=9,]
catch_area <- pfma %>% st_intersection(plot_area_ncc) %>% st_difference(bc_bound_hres())

bc_fish <- bc_bound_hres() %>% st_intersection(catch_area)
bc_map <- ggplot() +
  geom_sf(data = pfma, aes(fill=LABEL)) +
  geom_sf(data = bc_fish, fill = "grey10") +
  labs(fill="PFMA sub-areas") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(legend.position = "right")
bc_map


ncc_bec <- ggplot() +
  geom_sf(data = catch_area, aes(fill=as.factor(MGMT_AREA))) +
  scale_fill_manual(name="Management areas",values=c("#7570b3","#d95f02","#e7298a","#a6761d")) +
  new_scale("fill") +
  geom_sf(data=bec, aes(fill=ZONE,col=ZONE), alpha=0.8) +
  geom_sf(data = rivers_in_plot_area, colour = "lightblue1") +  #Plot Rivers
  geom_sf(data = rivers_ncc, colour = "lightblue1") +
  geom_sf(data = lakes_in_plot_area, fill = "lightblue1",colour=NA) +     #Plot Lakes
  geom_sf(data = coast_line, fill = NA) +                 #Plot coastline
  #scale_fill_viridis(name = "Last logging year",option="viridis",direction=-1) +
  #scale_fill_gradient2(name="Last logged year",midpoint = 1920, low="darkgreen",mid="olivedrab3",high="goldenrod1",space="Lab") +
  scale_fill_brewer(name="Biogeoclimatic zone",palette = "YlGn",direction=-1) +
  scale_colour_brewer(name="Biogeoclimatic zone",palette = "YlGn",direction=-1) +
  #scale_fill_brewer(name="Last logged year",palette = "Paired",direction=-1) +
  #ggsflabel::geom_sf_label_repel(data=data_point_labels[data_point_labels$population %in% c("Neekas River","Bella Coola River"),],aes(label=population),size=1.5, force = 1, nudge_x = -2, seed = 10)+ #Add labels
  ggsflabel::geom_sf_label_repel(data=data_point_labels,aes(label=population),size=1.5, force = 1, nudge_x = -2, seed = 10)+ #Add labels
  geom_sf(data=data_point_labels,pch=21,fill="white") +
  geom_sf(data = plot_area_ncc, alpha = 0,colour='black') +        #Plot area box
  coord_sf(expand = FALSE) +                                    #Expands box to axes
  xlab('Longitude') + ylab('Latitude') +                        #Axis titles
  annotation_scale(location = "br", width_hint = 0.5) +         #Rose Compass
  annotation_north_arrow(location = "br", which_north = "true",
                         pad_x = unit(0.5, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(1,"cm"), width = unit(1, "cm"))+
  theme(panel.background = element_rect('lightblue1'), panel.grid.major = element_line('lightblue1'),legend.position="top",legend.box.just="center",legend.box="horizontal",legend.justification = "center",legend.key.size=unit(1, "lines"),legend.margin = margin(c(0,0,0,0),unit="lines"),legend.title=element_text(size=7),legend.text = element_text(size=5))

bc_map <- ggplot() +
  geom_sf(data = bc_bound(), fill = "grey10") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(legend.justification = c(0, 1),legend.position = c(0, .95)) +
  theme(text = element_text(family = "Futura-Medium"),legend.title = element_text(family = "Futura-Bold", size = 10),legend.text = element_text(family = "Futura-Medium", size = 10))
bc <- bc_map + geom_rect(aes(xmin=ncc_extent@xmin,ymin=ncc_extent@ymin,xmax=ncc_extent@xmax,ymax=ncc_extent@ymax),fill = "tomato", colour = "grey70",alpha=0.5)

library(cowplot)
ncc_inset <- ggdraw(ncc_map) +
  draw_plot({
    bc},
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.8, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.70,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.15, 
    height = 0.15)

#Save the plot
ggsave('Figures/ncc_map.pdf',plot=ncc_inset,width = 6, height = 6,units='in',dpi=800)
ggsave('Figures/ncc_map.jpeg',plot=ncc_inset,width = 6, height = 6,units='in',dpi=800)
ncc_inset <- ggdraw(ncc_bec) +
  draw_plot({
    bc},
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.8, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.70,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.15, 
    height = 0.15)
ggsave('Figures/ncc_bec.pdf',plot=ncc_inset,width = 6, height = 6,units='in',dpi=800)
ggsave('Figures/ncc_bec.jpeg',plot=ncc_inset,width = 6, height = 6,units='in',dpi=800)
