

#################################    POSTERIOR PREDICTIVE CHECKS  #########################################


#SEED FOR REPRODUCTIBILITY
set.seed(42)


library("smoothSDE")
library("ggplot2")

#####################              CHECK GOODNESS OF FIT OF RESPONSE MODELS           ################################

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


check_response1=response1$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response2=response2$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response3=response3$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response4=response4$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response5=response5$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response6=response6$check_post(check_fn, n_sims = 500, silent = TRUE)



check_response1offset=response1offset$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response2offset=response2offset$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response3offset=response3offset$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response4offset=response4offset$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response5offset=response5offset$check_post(check_fn, n_sims = 500, silent = TRUE)
check_response6offset=response6offset$check_post(check_fn, n_sims = 500, silent = TRUE)






#######################   SIMULATE DATA WITH FITTED PARAMETERS MODEL AND  TRY TO RECOVER IT #########################




