
#' Code to process the geojson and shp files for greenland geography. 
#' We obtain a spatial features object (R package sf) that is more convenient to manipulate.



#to get polygons for the shore
#library(jsonlite)
library(sf)
library(ggplot2)
#library(ggmap)


#to deal with dataframes
library(dplyr) 


# Set the path to the directory containing coastline data
par_dir=dirname(dirname(getwd())) #parent directory 
data_path <- file.path(par_dir,"Data", "Greenland")

#get the coastline geometry from the geojson file
coastline<-st_read(file.path(data_path,"coastline.geojson"))
land<-st_read(file.path(data_path,"land.shp"))



#convert long lat coordinates to UTM zone 26
coastline_utm <- st_transform(coastline, crs = "+init=EPSG:32626 +units=km")
land_utm <- st_transform(land, crs = "+init=EPSG:32626 +units=km")


#crop again to actually get a rectangle in UTM coordinates
x.min= 410
x.max = 620
y.min = 7760
y.max = 7950

bbox <- st_bbox(c(xmin=x.min, xmax = x.max, ymin =y.min, ymax = y.max), 
                crs = st_crs(coastline_utm))

bbox_polygon <- st_as_sfc(st_bbox(bbox)
                          , crs = st_crs(coastline_utm))

land_utm<- land_utm %>% st_intersection(bbox_polygon)
test_utm<- test_utm %>% st_intersection(bbox_polygon)
coastline_utm<- coastline_utm %>% st_intersection(bbox_polygon)

#check polygons
ggplot()+geom_sf(data=coastline_utm$geometry,fill=NA)+
  coord_sf(datum=st_crs(32626))
ggplot()+geom_sf(data=land_utm$geometry,fill=NA)+
  coord_sf(datum=st_crs(32626))


#save the sf object in the Data directory
st_write(coastline_utm,file.path(data_path,"coastline_utm.geojson"),append=FALSE)
st_write(land_utm,file.path(data_path,"land_utm.shp"),append=FALSE)
