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

library(ggplot2)


# Read csv files ----------
files_DistanceShore<- list.files(pattern = "^result_rect_DistanceShore[0-9]+\\.csv$")
files_ExpShore<- list.files(pattern = "^result_rect_ExpShore[0-9]+\\.csv$")

df_ExpShore=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)
df_DistanceShore=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)

for (file in files_DistanceShore) {
  
  df=read.csv(file)
  df_DistanceShore=rbind(df_DistanceShore,df)
  
}

for (file in files_ExpShore) {
  
  df=read.csv(file)
  df_ExpShore=rbind(df_ExpShore,df)
  
}

# Bias and rmse for DistanceShore splines ------------

df_DistanceShore=df_DistanceShore[-1,]

N_DistanceShore=length(files_DistanceShore)

biais_DistanceShore=c()
re_biais_DistanceShore=c()
rmse_DistanceShore=c()
true_values_DistanceShore=c()

#loop over the coefficients
for (name in unique(df_DistanceShore$coeff_name)) {
  
  #all estimates for this coeff
  estimates=df_DistanceShore[df_DistanceShore$coeff_name==name,"estimate"]
  true_value=df_DistanceShore[df_DistanceShore$coeff_name==name,"true"][1]
  
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

results_DistanceShore=data.frame("coeff_name"=unique(df_DistanceShore$coeff_name),"true"=true_values_DistanceShore,
                                 "biais"=biais_DistanceShore,"re_biais"=re_biais_DistanceShore,"rmse"=rmse_DistanceShore)


p_DistanceShore=ggplot()+geom_violin(data=df_DistanceShore,aes(x=coeff_name,y=estimate,fill=coeff_name))+
  xlab(" ")+labs(fill = "Estimated coefficients")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_point(data = df_DistanceShore, aes(x = coeff_name, y = true), color = "red", size = 1,shape=4)

# Bias and rmse for DistanceShore splines ------------

res_ExpShore=res_ExpShore[-1,]

N_ExpShore=length(files_ExpShore)

biais_ExpShore=c()
re_biais_ExpShore=c()
rmse_ExpShore=c()
true_values_ExpShore=c()

#loop over the coefficients
for (name in unique(res_ExpShore$coeff_name)) {
  
  #all estimates for this coeff
  estimates=res_ExpShore[res_ExpShore$coeff_name==name,"estimate"]
  true_value=res_ExpShore[res_ExpShore$coeff_name==name,"true"][1]
  
  true_values_ExpShore=c(true_values_ExpShore,true_value)
  
  if (is.na(true_value)) {
    re_biais_ExpShore=c(re_biais_ExpShore,NA)
    biais_ExpShore=c(biais_ExpShore,NA)
    rmse_ExpShore=c(rmse_ExpShore,NA)
  }
  else {
    
    biais_ExpShore=c(biais_ExpShore,1/N_ExpShore*sum(estimates-true_value))
    rmse_ExpShore=c(rmse_ExpShore,sqrt(1/N_ExpShore*sum((estimates-true_value)^2)))
    
    if (true_value==0) {
      re_biais_ExpShore=c(re_biais_ExpShore,NA)
    }
    else {
      re_biais_ExpShore=c(re_biais_ExpShore,100*1/N_ExpShore*sum((estimates-true_value)/true_value))
    }
  }
}

results_ExpShore=data.frame("coeff_name"=unique(res_ExpShore$coeff_name),"true"=true_values_ExpShore,
                                 "biais"=biais_ExpShore,"re_biais"=re_biais_ExpShore,"rmse"=rmse_ExpShore)


p_ExpShore=ggplot()+geom_violin(data=df_ExpShore,aes(x=coeff_name,y=estimate,fill=coeff_name))+
  xlab(" ")+labs(fill = "Estimated coefficients")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_point(data = df_ExpShore, aes(x = coeff_name, y = true), color = "red", size = 1,shape=4)
