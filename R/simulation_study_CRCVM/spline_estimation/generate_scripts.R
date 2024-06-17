# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-03
#
# Script Name:    generate multiple scripts for estimations within three domains:
# rectangular, circular, fjords
#
# Script Description:
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

seeds=1:300
script_ExpShore<-readLines("base_script_ExpShore.R")
script_DistanceShore<-readLines("base_script_DistanceShore.R")

domains=c("rect","circ","fjords")

for (domain in domains) {
  
  #create directory to store the files for this domain
  if (!(dir.exists(domain))) {
    dir.create(domain)
  }
  
  # Find the line numbers containing "domain_name=", assuming there's only one
  line_ExpShore<- grep("domain_name=", script_ExpShore)
  line_DistanceShore <- grep("domain_name=", script_DistanceShore)
  
  # Change the next character
  script_ExpShore[line_ExpShore] <- paste('domain_name="',domain,'"',sep="")
  script_DistanceShore[line_DistanceShore] <-  paste('domain_name="',domain,'"',sep="")
  
  for (seed in seeds) {
    # Find the line numbers containing "seed=", assuming there's only one
    line_ExpShore<- grep("seed=", script_ExpShore)
    line_DistanceShore <- grep("seed=", script_DistanceShore)
    
    # Change the next character
    script_ExpShore[line_ExpShore] <- paste("seed=",seed)
    script_DistanceShore[line_DistanceShore] <-  paste("seed=",seed)
    
  
    # Rewrite the modified script back
    writeLines(script_ExpShore[-c(1:14)],file.path(domain,paste(domain,"_ExpShore",seed,".R",sep="")))
    writeLines(script_DistanceShore[-c(1:14)],file.path(domain,paste(domain,"_DistanceShore",seed,".R",sep="")))
  }
}

