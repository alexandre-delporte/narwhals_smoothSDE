# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2025 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2025-01-18
#
# Script Name:    
#
# Script Description: script to shrink land polygons of 100m to 
# get all observations in sea.
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(here)
library(sf)
library(ggplot2)

# Set the path to the directory containing the data
greenland_data_path <- here("Data","preprocessed_data","greenland")

#get the land geometry from shapefile
border<-st_read(file.path(greenland_data_path,"updated_scoresby_sound_utm.shp"))
border <- st_transform(border, crs = "+init=EPSG:32626 +units=km")

# Shrink polygons by 100 meters
land_shrunk <- st_buffer(border, dist = -0.2)

# Fix geometries if needed
land_shrunk <- st_make_valid(land_shrunk)

# Identify non-empty geometries
non_empty <- !st_is_empty(land_shrunk)

# Filter the sf object to keep only non-empty geometries
land_shrunk_clean <- land_shrunk[non_empty, ]

# Check the result
land_shrunk_clean

ggplot()+geom_sf(data=border$geometry,fill=NA,col="red")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_sf(data=land_shrunk$geometry,fill="grey")

#save the sf object in the preprocessed data directory
preprocessed_data_path <- here("Data","preprocessed_data","greenland")
st_write(land_shrunk_clean,file.path(preprocessed_data_path,"shrunk_updated_scoresby_sound_utm.shp"),append=FALSE)

