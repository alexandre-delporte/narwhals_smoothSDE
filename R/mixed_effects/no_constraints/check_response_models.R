

#################################    POSTERIOR PREDICTIVE CHECKS  #########################################


#SEED FOR REPRODUCTIBILITY
set.seed(42)


library("smoothSDE")
library("ggplot2")
source("/home/delporta/Documents/Recherche/Codes/smoothSDE/R/utility.R")

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
plot_checks(check_response1,"response1")

check_response2=response2$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response2,"response2")

check_response3=response3$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response3,"response3")

check_response4=response4$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response4,"response4")

check_response5=response5$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response5,"response5")

check_response6=response6$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response6,"response6")






check_response1offset=response1offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response1offset,"response1offset")

check_response2offset=response2offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response2offset,"response2offset")

check_response3offset=response3offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response3offset,"response3offset")

check_response4offset=response4offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response4offset,"response4offset")

check_response5offset=response5offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response5offset,"response5offset")

check_response6offset=response6offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response6offset,"response6offset")






#######################   SIMULATE DATA WITH FITTED PARAMETERS MODEL AND  TRY TO RECOVER IT #########################

# models without offset
check_fe_estimations(ctcrw=response2,model_name="response2",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response3,model_name="response3",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response4,model_name="response4",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)





#models with offset
check_fe_estimations(ctcrw=response2offset,model_name="response2offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response3offset,model_name="response3offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response4offset,model_name="response4offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)


check_fe_estimations(ctcrw=response8offset,model_name="response8offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response9offset,model_name="response9offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)




