# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-09-20
#
# Script Name:  process_boundaries.R
#
# Script Description:
# Code to process the shapefile downloaded from 
# https://osmdata.openstreetmap.de/data/land-polygons.html
# We only keep polygons for Scoresby Sound region, 
# and project it in UTM coordinate North zone 26

# SETUP ------------------------------------
cat("\014")                 
rm(list = ls())           



#to get root of git repo
library(here)

#to manipulate polygons
library(sf)

#plots
library(ggplot2)

#to deal with dataframes
library(dplyr) 



raw_data_path <- here("Data","raw_data","greenland")

land<-st_read(file.path(raw_data_path,"land_polygons.shp"))

#bounding box to crop the data to keep the region of interest in WGS84
bbox <- st_bbox(c(xmin=-30, xmax = -20, ymin =69, ymax = 72), 
                crs = st_crs(land))

# Crop the shapefile to the region of interest
land_crop <- st_crop(land, bbox)


#convert long lat coordinates to UTM zone 26
land_crop_utm <- st_transform(land_crop, crs = "+init=EPSG:32626 +units=km")

#crop again to actually get a rectangle in UTM coordinates
x.min= 410
x.max = 620
y.min = 7760
y.max = 7950

bbox <- st_bbox(c(xmin=x.min, xmax = x.max, ymin =y.min, ymax = y.max), 
                crs = st_crs(land_crop_utm))

bbox_polygon <- st_as_sfc(st_bbox(bbox)
                          , crs = st_crs(land_crop_utm))

land_crop_utm<- land_crop_utm %>% st_intersection(bbox_polygon)

#Visualize the polygons
plot(st_geometry(land_crop_utm))

#save the sf object in the preprocessed data directory
preprocessed_data_path <- here("Data","preprocessed_data","greenland")
st_write(land_crop_utm,file.path(preprocessed_data_path,"scoresby_sound_utm.shp"),append=FALSE)
