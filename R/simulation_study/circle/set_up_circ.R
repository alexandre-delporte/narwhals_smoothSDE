# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-07
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/simulation_study/circle/set_up_circ.R
#
# Script Description: set up for the simulation study in a circular domain
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space
set.seed(42)                #seed for reproducibility


library(ggplot2)            #plots
library(plotly)             #html plots
library(sf)                 #geometry for constraints
library(mgcv)               #tensor splines


N_ID=6                      #number of individual tracks per batch
TMAX=1*24                 #duration of each track in hour
SP_DF=3                     #degree of freedom in tensor splines


# Cicular domain --------------------

# Create circle geometries
large_circle <- st_buffer(st_sfc(st_point(c(0,0))), dist = 10)
small_circle <- st_buffer(st_sfc(st_point(c(0,0))), dist = 8)

# Convert circles to polygons
large_circle_polygon <- st_cast(large_circle, "POLYGON")
small_circle_polygon <- st_cast(small_circle, "POLYGON")

# Create sf objects
large_circle_sf <- st_sf(geometry = large_circle_polygon)
small_circle_sf <- st_sf(geometry = small_circle_polygon)

# Calculate the symmetric difference to get the border
border <- st_sym_difference(large_circle_sf, small_circle_sf)
