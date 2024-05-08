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

files_DistanceShore<- list.files(pattern = "^result_rect_DistanceShore[0-9]+\\.csv$")
files_ExpShore<- list.files(pattern = "^result_rect_ExpShore[0-9]+\\.csv$")

res_ExpShore=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)
res_DistanceShore=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)

for (file in files_DistanceShore) {
  
  df=read.csv(file)
  res_DistanceShore=rbind(res_DistanceShore,df)
  
}

for (file in files_ExpShore) {
  
  df=read.csv(file)
  res_ExpShore=rbind(res_ExpShore,df)
  
}

res_ExpShore=res_ExpShore[-1,]
res_DistanceShore=res_DistanceShore[-1,]

N_DistanceShore=length(files_DistanceShore)
N_ExpShore=length(files_ExpShore)

biais_DistanceShore=c()
re_biais_DistanceShore=c()
rmse_DistanceShore=c()
true_values_DistanceShore=c()

#loop over the coefficients
for (name in unique(res_DistanceShore$coeff_name)) {
  
  #all estimates for this coeff
  estimates=res_DistanceShore[res_DistanceShore$coeff_name==name,"estimate"]
  true_value=res_DistanceShore[res_DistanceShore$coeff_name==name,"true"][1]
  
  true_values_DistanceShore=c(true_values_DistanceShore,true_value)
  
  if (is.na(true_value)) {
    re_biais_DistanceShore=c(re_biais_DistanceShore,NA)
    biais_DistanceShore=c(biais_DistanceShore,NA)
    rmse_DistanceShore=c(rmse_DistanceShore,NA)
  }
  else {
    
    biais_DistanceShore=c(biais_DistanceShore,1/N_DistanceShore*sum(estimates-true_value))
    rmse_DistanceShore=c(rmse_DistanceShore,sqrt(1/N_DistanceShore*sum((estimates-true_value)^2)))
    
    if (true_value==0) {
      re_biais_DistanceShore=c(re_biais_DistanceShore,NA)
    }
    else {
      re_biais_DistanceShore=c(re_biais_DistanceShore,100*1/N_DistanceShore*sum((estimates-true_value)/true_value))
    }
  }
}
  
names(true_values_DistanceShore)=unique(res_DistanceShore$coeff_name)
names(biais_DistanceShore)=unique(res_DistanceShore$coeff_name)
names(re_biais_DistanceShore)=unique(res_DistanceShore$coeff_name)
names(rmse_DistanceShore)=unique(res_DistanceShore$coeff_name)




                      
