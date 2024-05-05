# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-05
#
# Script Name:    get_results.R
#
# Script Description: Compute biais and rmse from the results in the .csv files
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

res_ExpShore=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)
res_DistanceShore=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)
for (k in 1:100) {
  
  df_ExpShore=read.csv(paste("result_ExpShore",k,".R",sep=""))
  df_DistanceShore=read.csv(paste("result_DistanceShore",k,".R",sep=""))
  
  res_ExpShore=rbind(res_ExpShore,df_ExpShore)
  res_DistanceShore=rbind(res_DistanceShore,df_DistanceShore)
  
}
