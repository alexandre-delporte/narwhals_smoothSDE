# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-03
#
# Script Name:    
#
# Script Description:
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

seeds=1:100
script_ExpShore<-readLines("base_script_rect_ExpShore.R")
script_DistanceShore<-readLines("base_script_rect_DistanceShore.R")

for (seed in seeds) {
  # Find the line numbers containing "seed=", assuming ther's only one
  line_ExpShore<- grep("seed=", script_ExpShore)
  line_DistanceShore <- grep("seed=", script_DistanceShore)
  
  # Change the next character
  script_ExpShore[line_ExpShore] <- paste("seed=",seed)
  script_DistanceShore[line_ExpShore] <-  paste("seed=",seed)
  
  # Rewrite the modified script back
  writeLines(script_ExpShore[-c(1:13)], paste("rect_ExpShore",seed,".R",sep=""))
  writeLines(script_DistanceShore[-c(1:13)], paste("rect_DistanceShore",seed,".R",sep=""))
}

