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
# Script Description: script to shrink land polygons to 
# get all observations in sea.
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(here)
library(sf)
library(ggplot2)
library(plotly)
library(htmlwidgets)

# Set the path to the directory containing the data
greenland_data_path <- here("Data","preprocessed_data","greenland")

#get the land geometry from shapefile
border<-st_read(file.path(greenland_data_path,"updated_scoresby_sound_utm.shp"))
border <- st_transform(border, crs = "+init=EPSG:32626 +units=km")


shrink_border=function(border,dist) {
  # Shrink polygons by 100 meters
  land_shrunk <- st_buffer(border, dist = -dist)

  # Fix geometries if needed
  land_shrunk <- st_make_valid(land_shrunk)

  # Identify non-empty geometries
  non_empty <- !st_is_empty(land_shrunk)

  # Filter the sf object to keep only non-empty geometries
  land_shrunk_clean <- land_shrunk[non_empty, ]
  
  land_shrunk_clean <- st_cast(land_shrunk_clean, "MULTIPOLYGON")

  # Check the result
  land_shrunk_clean
}

land_shrunk100<-shrink_border(border,dist=0.1)
land_shrunk200<-shrink_border(border,dist=0.2)
land_shrunk300<-shrink_border(border,dist=0.3)


# Plot the three adjusted borders
land_shrunk100$distance <- "100m"
land_shrunk200$distance <- "200m"
land_shrunk300$distance <- "300m"
border$distance <- "0m"

land_shrunk_combined <- rbind(land_shrunk100, land_shrunk200, land_shrunk300,border)

# Plot with legend
plot_borders <- ggplot() +
  geom_sf(data = land_shrunk_combined, aes(col = distance), fill = "grey") +
  scale_color_manual(values = c("0m"="red","100m" = "orange", "200m" = "green", "300m" = "blue"),
                    name = "Distance from Shore") + 
  coord_sf(crs = 32626) +
  theme_minimal()

dyn_plot_borders<-ggplotly(plot_borders)

saveWidget(dyn_plot_borders,"shrunk_borders.html")

#save the sf object in the preprocessed data directory
preprocessed_data_path <- here("Data","preprocessed_data","greenland")
st_write(land_shrunk100,file.path(preprocessed_data_path,
                                  "shrunk100_scoresby_sound_utm.shp"),append=FALSE)
st_write(land_shrunk200,file.path(preprocessed_data_path,
                                  "shrunk200_scoresby_sound_utm.shp"),append=FALSE)
st_write(land_shrunk300,file.path(preprocessed_data_path,
                                  "shrunk300_scoresby_sound_utm.shp"),append=FALSE)

