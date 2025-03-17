# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-03
#
# Script Name:    generate scripts for estimations within three domains:
# rectangular, circular, fjords with hyperparameters set up defined in the scripts 
# 'set_uo_rect','set_up_circ','set_up_fjords'
#
# Script Description:
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

N_SCRIPTS=100
seeds=1:N_SCRIPTS
base_script=paste0("base_script.R")
script<-readLines(base_script)

domain="fjords"


  
#create directory to store the files for this domain
if (!(dir.exists(file.path(domain,"Rscripts")))) {
  dir.create(file.path(domain,"Rscripts"))
}

# Find the line numbers containing "domain_name=", assuming there's only one
line<- grep("domain_name=", script)
  
# Change the next character
script[line] <- paste('domain_name="',domain,'"',sep="")
  
  for (seed in seeds) {
    # Find the line numbers containing "seed=", assuming there's only one
    line<- grep("seed=", script)
    
    # Change the next character
    script[line] <- paste("seed=",seed)
    
  
    # Rewrite the modified script back
    writeLines(script[-c(1:14)],file.path(domain,"Rscripts",paste(domain,seed,".R",sep="")))
  }


