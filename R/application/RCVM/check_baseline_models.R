
# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-30
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/application/CRCVM/check_baseline_models.R
#
# Script Description:
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console

#SEED FOR REPRODUCTIBILITY
set.seed(42)


#####################              CHECK GOODNESS OF FIT OF BASELINE MODEL            ################################

#' @description compute mean of empirical velocity phase and norm for each individual track
#' @return dataframe with one column per individual and one row per statistic 
check_fn=function(data) {
  
  n_id=length(unique(data$ID))
  stats_df=as.data.frame(matrix(rep(NA,n_id*2),ncol=n_id))
  colnames(stats_df)=unique(data$ID)
  rownames(stats_df)=c("mean_vnorm","mean_phase")
  
  for (id in unique(data$ID)) {
    
    #data for specific id
    sub_data=data[data$ID==id,]
    n_sub=length(sub_data[,1])
    
    #time steps
    dtimes=sub_data[2:n_sub,"time"]-sub_data[1:(n_sub-1),"time"]
    
    #step lengths
    dx=(sub_data[2:n_sub,"x"]-sub_data[1:(n_sub-1),"x"])
    dy=(sub_data[2:n_sub,"y"]-sub_data[1:(n_sub-1),"y"])
    
    #empirical velocity norms
    vexp_df=cbind(dx/dtimes,dy/dtimes)
    vnorm=sqrt(vexp_df[,1]^2+ vexp_df[,2]^2)
    
    #empirical phase angles
    phase=signed_angle(matrix(rep(c(1,0),each=n_sub-1),ncol=2),as.matrix(vexp_df))
    
    stats_df[,id]=c(mean(vnorm),mean(phase))
  }
  
  return (stats_df)
}

check_baseline3=baseline3$check_post(check_fn, model_name="baseline3",n_sims =500, silent = TRUE)



