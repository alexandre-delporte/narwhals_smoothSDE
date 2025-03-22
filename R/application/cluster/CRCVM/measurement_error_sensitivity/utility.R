

#' Plot interpolated trajectories, angles and distances
#' @param data Dataframe with the columns time, ID, BoundaryDistance, and BoundaryAngle
#' @param interpolation_data 
#' @param border Boundary of the domain as a sf object with geometry 
#' @param response name of response variables
#' @return List of plots
get_interpolation_plots<-function(data,interpolation_data,border,response=c("x","y")) {
  
  
  #create dataframe with interpolated values
  end_indices <- c(which(diff(as.numeric(factor(data$ID))) != 0),nrow(data))
  n_step<-ncol(interpolation_data$BoundaryDistance)
  interpolated_x_vec=as.vector(t(interpolation_data[[response[1]]][-end_indices,]))
  interpolated_y_vec=as.vector(t(interpolation_data[[response[2]]][-end_indices,]))
  interpolated_distance_vec=as.vector(t(interpolation_data$BoundaryDistance[-end_indices,]))
  interpolated_angle_vec=as.vector(t(interpolation_data$BoundaryAngle[-end_indices,]))
  interpolation_time=as.vector(t(interpolation_data$time[-end_indices,]))
  ID=unlist(lapply(unique(data$ID), function(id) {
    n=nrow(data[data$ID==id,])
    rep(id,(n-1)*n_step)
  }))
  
  #combine dataframes
  interpolated_df=data.frame("time"=interpolation_time,"ID"=ID,"x"=interpolated_x_vec,
                             "y"=interpolated_y_vec,
                             "BoundaryDistance"=interpolated_distance_vec,
                             "BoundaryAngle"=interpolated_angle_vec)
  interpolated_df$type="Interpolated"
  data$type="Observed"
  colnames(data)[colnames(data)==response[1]]<-"x"
  colnames(data)[colnames(data)==response[2]]<-"y"
  
  combined_data<-rbind(interpolated_df,data)
  
  
  # Create plot with legends for shape and linetype
  trajectory_plot = ggplot() +geom_sf(data=border$geometry,fill="grey")+
    geom_point(data = combined_data, 
               aes(x,y, col = ID, shape = type), size = 0.2,alpha=1) +
    
    # Customize the legend
    labs(color = "Trajectory ID", 
         shape = "Point Type", 
         linetype = "Line Type", 
         title = "Interpolated vs. Observed trajectories") +
    scale_shape_manual(values = c(10:12))+
    theme_minimal() +
    theme(legend.position = "right")  # Adjust legend position
  
  
  
  distance_plot = ggplot() +
    # Scatter plot for observed BoundaryDistance (Shape Legend)
    geom_line(data = combined_data, 
              aes(time, BoundaryDistance, col = ID, linetype = type), size = 0.5) +
    
    # Customize the legend
    labs(color = "Trajectory ID", 
         shape = "Point Type", 
         linetype = "Line Type", 
         x = "Time", y = "Boundary Distance", 
         title = "Interpolated vs. True Boundary Distance") +
    
    theme_minimal() +
    theme(legend.position = "right")  # Adjust legend position
  
  
  # Create plot with legends for shape and linetype
  angle_plot = ggplot() +
    # Scatter plot for observed BoundaryDistance (Shape Legend)
    geom_line(data = combined_data, 
              aes(time, BoundaryAngle, col = ID, linetype = type), size = 0.5) +
    
    # Customize the legend
    labs(color = "Trajectory ID", 
         shape = "Point Type", 
         linetype = "Line Type", 
         x = "Time", y = "Boundary Angle", 
         title = "Interpolated vs. True Boundary Angle") +
    
    theme_minimal() +
    theme(legend.position = "right")  # Adjust legend position
  
  return(list("angle_plot"=angle_plot,"distance_plot"=distance_plot,"trajectory_plot"=trajectory_plot))
}


#' Fit baseline on data
#' @param data data before exposure 
#' @param interpolation_data 
#' @param sigma_obs GPS measurement error std
#' @param a,b,D0,D1,sigma_D,sigma_theta parameters in the function omega
#' @param fixpar parameters to fix
#' @param response name of response variables
#' @return Fitted baseline SDE object
fit_baseline<-function(data,interpolation_data,sigma_obs=0.05,a=2,b=1,D0=0.3,D1=0.5,
                       sigma_D=1,sigma_theta=0.6,fixpar=c("a","b","D0","D1","sigma_D","sigma_theta"),
                       response=c("x","y")) {
  
  #Fixed measurement error
  H=array(rep(sigma_obs^2*diag(2),length(data$time)),dim=c(2,2,length(data$time)))
  
  interpolated_dist<- interpolation_data$BoundaryDistance
  interpolated_dist[interpolated_dist == 0] <- 0.1
  
  #parameters formulas
  formulas <- list(tau =~s(ID,bs="re"),nu=~s(ID,bs="re"),
                   a=~1,b=~1,D0=~1,D1=~1,sigma_D=~1,sigma_theta=~1)
  
  # Fit baseline
  baseline <- SDE$new(formulas = formulas,
                      data = data,type = "CRCVM_SSM",
                      response = response,par0 = c(1.2,4.5,a,b,D0,D1,sigma_D,sigma_theta),
                      fixpar=fixpar,
                      other_data=
                        list("H"=H,"interpolated_distance"=interpolated_dist,
                             "interpolated_angle"=interpolation_data$BoundaryAngle))
  baseline$update_lambda(c(1,1))
  baseline$fit()
  
  return (baseline)
}



#' Get estimated from fitted baseline SDE model
#' @param sde fitted sde
#' @param include_sigma_obs boolean
#' @param n_post number of draws from the posterior distribution used to give estimates
#' and confidence intervals
#' @param resp Logical (default: TRUE). Should the output be on 
#' the response scale? If FALSE, the output is on the linear 
#' predictor scale
#' @return List of estimates, standard error and quantiles for each parameter
get_baseline_estimates<- function(sde,include_sigma_obs=FALSE,n_post=100,
                                  resp=FALSE) {
  
  post_coeff=sde$post_coeff(n_post=n_post)
  link=ifelse(resp,exp,identity)
  post_par=list(
    "tau"=link(post_coeff$coeff_fe[,"tau.(Intercept)"]),
    "nu"=link(post_coeff$coeff_fe[,"nu.(Intercept)"]),
    "sigma_tau"=1/sqrt(exp(post_coeff$log_lambda[,"tau.s(ID)"])),
    "sigma_nu"=1/sqrt(exp(post_coeff$log_lambda[,"nu.s(ID)"])))
  
  if (include_sigma_obs) {
    post_par$sigma_obs=1000*exp(post_coeff$log_sigma_obs)
  }
  
  mean_par=lapply(post_par,mean)
  sd_par=lapply(post_par,sd)
  quant_par <- lapply(post_par,
                      quantile, probs = c(0.025, 0.975))
  
  return(list("estimates"=mean_par,"std.error"=sd_par,"quantiles"=quant_par))
  
}

#' Make latex table to show the estimates and confidence intervals
#' @param sde fitted sde
#' @param include_sigma_obs boolean
#' @param n_post number of draws from the posterior distribution used to give estimates
#' and confidence intervals
#' @return NULL (just print latex table)
make_final_baseline_tab<-function(sde,include_sigma_obs=FALSE,n_post=100,resp=FALSE) {
  
  est=get_baseline_estimates(sde,include_sigma_obs,n_post,resp)
  mean_par=est$estimates
  sd_par=est$std.error
  quant_par=est$quantiles
  parameter_names <- names(mean_par)
  estimates <- mean_par
  conf_intervals <- sapply(parameter_names, function(param) {
    sprintf("$[%.2f; %.2f]$", quant_par[[param]][1], 
            quant_par[[param]][2])
  })
  
  parameter_names=c("$\\tau$", "$\\nu$","$\\sigma_{\\tau}$","$\\sigma_{\\nu}$")
  if (include_sigma_obs) {
    parameter_names=c(parameter_names,"$\\sigma_{\\obs}$")
  }
  table_data <- data.frame(
    Parameter = parameter_names,
    Estimate = sprintf("$%.2f$", mean_par),
    CI = conf_intervals,
    stringsAsFactors = FALSE)
  
  xtab <- xtable(
    table_data,
    caption = paste(sde$type(),"baseline estimations with random effects"),
    label = paste0("table:",sde$type(),"_baseline_estimations_with_random_effects"))
  
  print(
    xtab,type="latex",
    include.rownames = FALSE,
    sanitize.text.function = identity,  # Prevent escaping of LaTeX math symbols
    hline.after = c(-1, 0, nrow(table_data)),  
    table.placement="H")
  
} 

#' Make plot to compare distributions in simulation and observation
#' @param simulated_data list of data frames simulated from an SDE model, each containing columns "BoundaryDistance", "BoundaryAngle" and "speed"
#' @param observed_data real GPS data (data frame).
#' @return ggplot of boundary distances distributions
make_density_plots <- function(simulated_data, observed_data) {
  
  # Ensure simulated_data is a list of data frames
  if (!is.list(simulated_data) || !all(sapply(simulated_data, is.data.frame))) {
    stop("simulated_data must be a list of data frames")
  }
  
  
  simulated_combined <- bind_rows(
    lapply(seq_along(simulated_data), function(i) {
      simulated_data[[i]] %>% mutate(SimulationID = as.factor(i), Type = "Simulated")
    })
  )
  
  observed_combined <- observed_data %>% mutate(Type = "Observed")
  
  combined_data <- bind_rows(simulated_combined, observed_combined)
  
  # Plot density distribution
  plot_distance <- ggplot(combined_data, aes(x = BoundaryDistance, color = Type, fill = Type)) +
    geom_density(data = simulated_combined, aes(group = SimulationID), alpha = 0.1) + 
    geom_density(data = observed_combined, alpha = 0.2) +
    xlab("Distance to Shore") +
    scale_color_manual(values = c("Simulated" = "blue", "Observed" = "red")) + 
    scale_fill_manual(values = c("Simulated" = "blue", "Observed" = "red")) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.title.x = element_blank())+
  theme(legend.title = element_blank())
  
  
  plot_angle <- ggplot(combined_data, aes(x = abs(BoundaryAngle), color = Type, 
                                          fill = Type)) +
    geom_density(data = simulated_combined, aes(group = SimulationID), 
                 alpha = 0.1) + 
    geom_density(data = observed_combined, alpha = 0.2) + 
    xlab(expression(Theta)) +
    scale_color_manual(values = c("Simulated" = "blue", "Observed" = "red")) + 
    scale_fill_manual(values = c("Simulated" = "blue", "Observed" = "red")) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.title.x = element_blank())+
  theme(legend.title = element_blank())
  
  
  plot_speed <- ggplot(combined_data, aes(x = speed, color = Type, 
                                          fill = Type)) +
    geom_density(data = simulated_combined, aes(group = SimulationID), 
                 alpha = 0.1) + 
    geom_density(data = observed_combined, alpha = 0.2) + 
    xlab("Speed") +
    scale_color_manual(values = c("Simulated" = "blue", "Observed" = "red")) + 
    scale_fill_manual(values = c("Simulated" = "blue", "Observed" = "red")) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.title.x = element_blank())+
  theme(legend.title = element_blank())
  
  
  
  return(list("distance"=plot_distance,"speed"=plot_speed,"angle"=plot_angle))
}

#' Compute empirical speed from finite differences
#' @param data data with column time
#' @param response name of response variables
#' @return vector of empirical speeds
compute_speed<-function(data,response=c("x","y")) {
  
  speed_vector<-rep(0,nrow(data))
  
  for (id in unique(data$ID)){
    
    #filter data with specific id
    sub_ind=(data$ID==id)
    sub_data=data[sub_ind,]
    n_sub=length(sub_data$time)
    
    #time steps
    dtimes=sub_data[2:n_sub,"time"]-sub_data[1:(n_sub-1),"time"]
    
    #step lengths
    dx=sub_data[2:n_sub,response[1]]-sub_data[1:(n_sub-1),response[1]]
    dy=sub_data[2:n_sub,response[2]]-sub_data[1:(n_sub-1),response[2]]
    
    speed=c(NA,sqrt((dx/dtimes)^2+(dy/dtimes)^2))
    
    speed_vector[sub_ind]<-speed
  }
  
  return(speed_vector)
}

plot_trajectories=function(simulated_trajectories,observed_trajectories) {
  
  plot_list=list()
  
  
  for (i in seq_along(simulated_trajectories)) {
  
    plot_sim_trajectory=ggplot() +
      geom_sf(data=border200$geometry,fill="grey")+
      coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
      geom_path(data=simulated_trajectories[[i]], 
                aes(x = x, y = y,color=ID),size = 0.2, 
                lineend = "round",alpha=0.5) +     
      geom_point(data=simulated_trajectories[[i]],
                 aes(x = x, y = y,color=ID),size =0.2 ) +
      theme_minimal()+
      scale_color_viridis_d(name="ID",labels=c("A1","A2","A3","A4","A5","A6"))+theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())
    
    plot_obs_trajectory=ggplot() +
      geom_sf(data=border200$geometry,fill="grey")+
      coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
      geom_path(data=observed_trajectories, 
                aes(x = x, y = y,color=ID),size = 0.2, 
                lineend = "round",alpha=0.5) +     
      geom_point(data=observed_trajectories,
                 aes(x = x, y = y,color=ID),size =0.2 ) +
      theme_minimal()+
      scale_color_viridis_d(name="ID",labels=c("A1","A2","A3","A4","A5","A6"))+theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())
    
    
    final_plot=ggarrange(plotlist=list(plot_obs_trajectory,plot_sim_trajectory),
                         ncol=2, nrow=1, common.legend = TRUE, legend="bottom",
                         labels= c("Observed","Simulated"), vjust=2)
    
    plot_list[[i]]=final_plot
  }
  
  return (plot_list)
}



#' Fit response model on data
#' @param baseline fitted sde model before exposure
#' @param data data after exposure 
#' @param interpolation_data 
#' @param sigma_obs GPS measurement error std
#' @return Fitted baseline SDE object
fit_response<-function(baseline,data,interpolation_data,sigma_obs=0.05)  {
  
  
  #Fixed measurement error
  H=array(rep(sigma_obs^2*diag(2),nrow(data)),dim=c(2,2,nrow(data)))
  
  interpolated_distance<- interpolation_data$BoundaryDistance
  interpolated_distance[interpolated_distance == 0] <- 0.1
  
  #parameters formulas
  formulas <- list(tau =~ExpShip+s(ID,bs="re"),nu=~ExpShip+s(ID,bs="re"),
                   a=~1,b=~1,D0=~1,D1=~1,sigma_D=~1,sigma_theta=~1)
  
  # Create SDE object
  response_crcvm<- SDE$new(formulas = formulas,
                           data = data,type = "CRCVM_SSM",
                           response = baseline$response(),par0 = baseline$par(),
                           fixpar=c("a","b","D0","D1","sigma_theta","sigma_D"),
                           other_data=
                             list("H"=H,"interpolated_distance"=interpolated_distance,
                                  "interpolated_angle"=interpolation_data$BoundaryAngle))
  
  
  new_map=list(coeff_re=factor(rep(NA,12)),coeff_fe=factor(c(NA,1,NA,2,rep(NA,6))),
               log_lambda=factor(c(NA,NA)))
  
  new_coeff_fe=c(baseline$coeff_fe()["tau.(Intercept)",1],0,baseline$coeff_fe()["nu.(Intercept)",1],0,
                 baseline$coeff_fe()[3:8,1])
  new_coeff_re=baseline$coeff_re()
  new_lambda=c(baseline$lambda())
  
  response_crcvm$update_map(new_map)
  response_crcvm$update_coeff_re(new_coeff_re)
  response_crcvm$update_coeff_fe(new_coeff_fe)
  response_crcvm$update_lambda(new_lambda)
  
  response_crcvm$fit()
  
  return (response_crcvm)
}


#' Fit response model using prosterior draws of the baseline coefficients
#' @param baseline fitted sde model before exposure
#' @param data data after exposure 
#' @param interpolation_data
#' @param n_sim number of posterior draws from the baseline
#' @param sigma_obs GPS measurement error
propagate_baseline_uncertainty<-function(baseline,data,interpolation_data,n_sim=100,sigma_obs=0.05) {
  
  #posterior samples of the coefficient of the baseline model
  post_coeff_bas=baseline$post_coeff(n_post=n_sim)
  post_coeff_re=post_coeff_bas$coeff_re
  post_coeff_fe=post_coeff_bas$coeff_fe
  post_coeff_lambda=exp(post_coeff_bas$log_lambda)
  
  #define model formulas
  formulas <- list(tau =~ExpShip+s(ID,bs="re"),
                   nu=~ExpShip+s(ID,bs="re"),a=~1,b=~1,D0=~1,
                   D1=~1,sigma_D=~1,sigma_theta=~1)
  
  #Fixed measurement error
  H=array(rep(sigma_obs^2*diag(2),nrow(data)),dim=c(2,2,nrow(data)))
  
  interpolated_distance<- interpolation_data$BoundaryDistance
  interpolated_distance[interpolated_distance == 0] <- 0.1
  
  # vector to store estimates and std
  alpha_estimates=matrix(rep(0,2*nrow(post_coeff_fe)),ncol=2)
  alpha_std=matrix(rep(0,2*nrow(post_coeff_fe)),ncol=2)
  
  #list to store models
  models_list=list()
  
  for (i in 1:nrow(post_coeff_fe)) {
    
    
    
    response_crcvm_temp<- SDE$new(formulas = formulas,
                                  data = data,type = "CRCVM_SSM",
                                  response = baseline$response(),par0 = baseline$par(),
                                  fixpar=c("a","b","D0","D1","sigma_theta","sigma_D"),
                                  other_data=
                                    list("H"=H,"interpolated_distance"=interpolated_distance,
                                         "interpolated_angle"=interpolation_data$BoundaryAngle))
    
    
    new_map=list(coeff_re=factor(rep(NA,12)),coeff_fe=factor(c(NA,1,NA,2,rep(NA,6))),
                 log_lambda=factor(c(NA,NA)))
    
    new_coeff_fe=c(post_coeff_fe[i,"tau.(Intercept)"],0,post_coeff_fe[i,"nu.(Intercept)"],0,
                   post_coeff_fe[i,3:8])
    new_coeff_re=post_coeff_re[i,]
    new_lambda=post_coeff_lambda[i,]
    
    response_crcvm_temp$update_map(new_map)
    response_crcvm_temp$update_coeff_re(new_coeff_re)
    response_crcvm_temp$update_coeff_fe(new_coeff_fe)
    response_crcvm_temp$update_lambda(new_lambda)
    
    try(response_crcvm_temp$fit())
    
    est_alpha=as.list(response_crcvm_temp$tmb_rep(),what="Est")$coeff_fe[c("tau.ExpShip","nu.ExpShip"),1]
    std_alpha=as.list(response_crcvm_temp$tmb_rep(),what="Std")$coeff_fe[c("tau.ExpShip","nu.ExpShip"),1]
    
    alpha_estimates[i,]=est_alpha
    alpha_std[i,]=std_alpha
    
    models_list[[i]]=response_crcvm_temp
  }
  
  colnames(alpha_estimates)=c("tau.ExpShip","nu.ExpShip")
  colnames(alpha_std)=c("tau.ExpShip","nu.ExpShip")
  
  return(list("models"=models_list,"estimates"=alpha_estimates,"std"=alpha_std))
}

