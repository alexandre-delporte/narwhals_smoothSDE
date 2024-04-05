
# NECESSARY PACKAGES

#for sde fit
library(smoothSDE)

#for plots
library(ggplot2)
library(plotly)

#to define and extract string formulas
library(splines)
library(formula.tools)
library(stringr)

#for gam models
library(mgcv)

#to save 3D files
library(htmlwidgets)

#for nearest point on coastline
library(sf)

#for block matrix in SSM for simulation
library(Matrix)
library(MASS)


####################### FUNCTION TO GET BEST INITIAL PARAMETERS ##########################


get_best_par0=function(sde,par0_df,fixpar=NULL) {
  
  formulas <- sde$formulas()
  data=sde$data()
  response=sde$response()
  type=sde$type()
  
  n=length(par0_df[,1])
  
  minAIC=Inf
  best_par0=c(0,0,0,0)
  
  for (i in 1:n) {
    par0 <- as.numeric(par0_df[i,])
    print(par0)
    baseline <- SDE$new(formulas = formulas,data = data,
                        type = type,response = response,par0 = par0,fixpar=fixpar)
    baseline$fit()
    
    AIC=baseline$AIC_marginal()
    #print(AIC)
    if (!is.nan(AIC) & AIC < minAIC) {
      best_par0=par0
      minAIC=AIC
    }
  }
  return (best_par0)
  
}

######################################################## FUNCTIONS TO CHECK MODELS ESTIMATIONS #############################################################################

#' @description compute mean of empirical velocity phase and norm for each individual track
#' @return dataframe with n_id columns and 2 rows 
check_fn=function(data) {
  
  n_id=length(unique(data$ID))
  stats_df=data.frame(rep(NA,n_id*2),ncol=n_id)
  
  for (id in unique(data$ID)) {
    
    #data for specific id
    sub_data=data[data$ID==id,]
    n=length(sub_data[,1])
    
    #time steps
    dtimes=sub_data[2:n,"time"]-sub_data[1:(n-1),"time"]
    
    #step lengths
    dx=(sub_data[2:n,"x"]-sub_data[1:(n-1),"x"])
    dy=(sub_data[2:n,"y"]-sub_data[1:(n-1),"y"])
    
    #empirical velocity norms
    vexp_df=cbind(dx/dtimes,dy/dtimes)
    vnorm=sqrt(vexp_df[,1]^2+ vexp_df[,2]^2)
    
    #empirical phase angles
    phase=signed_angle(matrix(rep(c(1,0),each=n_sub-1),ncol=2),vexp_df)
      
    stats_df[,id]=c(mean(vnorm),mean(phase))
    }
    
  return (stats_df)
}


plot_checks=function(check_stats,model_name) {
  
  n_id=6
  n_stats=3
  
  stats=check_stats$stats
  
  obs_stat=check_stats$obs_stat
  
  # Data frames for ggplot (observed and simulated values)
  df_obs <- as.data.frame.table(as.matrix(obs_stat))
  colnames(df_obs) <- c("stat", "sim", "val")
  df <- as.data.frame.table(stats)
  colnames(df) <- c("stat", "sim", "val")
  
  index=seq(from=1,to=length(df[,1]),by=n_stats)
  df_mean=df[index,]
  df_sd=df[index+1,]
  df_theta=df[index+2,]
  df_obs_mean=df_obs[seq(from=1,to=n_id*n_stats,by=n_stats),]
  df_obs_sd=df_obs[seq(from=1,to=n_id*n_stats,by=n_stats)+1,]
  df_obs_theta=df_obs[seq(from=1,to=n_id*n_stats,by=n_stats)+2,]
  
  # Create plot
  names=c(rep(paste("A",1:n_id,sep=""),each=n_stats))
  names(names)=paste("statistic",1:(n_id*n_stats),sep=" ")
  
  pmean <- ggplot(df_mean, aes(val)) + 
    geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                   bg = "lightgrey", na.rm = TRUE) +
    facet_wrap(~stat, scales = "free",labeller=as_labeller(names)) + 
    geom_vline(aes(xintercept = val), data = df_obs_mean,col="red") +
    theme_light() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black"))
  
  psd <- ggplot(df_sd, aes(val)) + 
    geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                   bg = "lightgrey", na.rm = TRUE) +
    facet_wrap(~stat, scales = "free",labeller=as_labeller(names)) + 
    geom_vline(aes(xintercept = val), data = df_obs_sd,col="red") +
    theme_light() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black"))
  
  ptheta <- ggplot(df_theta, aes(val)) + 
    geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                   bg = "lightgrey", na.rm = TRUE) +
    facet_wrap(~stat, scales = "free",labeller=as_labeller(names)) + 
    geom_vline(aes(xintercept = val), data = df_obs_theta,col="red") +
    theme_light() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black"))
  
  
  if (!(dir.exists(model_name))) {
    dir.create(model_name)
  }
  
  ggsave(paste(paste("mean_velocity_check",model_name,sep="_"),".png",sep=""),plot=pmean,path=model_name)
  ggsave(paste(paste("sd_velocity_check",model_name,sep="_"),".png",sep=""),plot=psd,path=model_name)
  ggsave(paste(paste("mean_angle_check",model_name,sep="_"),".png",sep=""),plot=ptheta,path=model_name)
}




check_fe_estimations=function(sde,model_name,links=list(),xmins=list(),xmaxs=list(),xlabels=list(),npost=1000,level=0.95) {
  
  #simulate a track without measurement errors from a fitted model
  #then fit the same model on this track 
  #and plot the estimated fixed effects with confidence interval vs the actual population smooth
  
  
  #get data
  data=sde$data()
  
  # get formulas and parameters
  formulas=sde$formulas()
  pars=names(formulas)
  n_par=length(formulas)
  
  #fixed coeffs (offsets)
  map=sde$map()
  
  #initial value of parameters 
  par0=sde$par()
  
  #simulate from the fitted model
  sde_sim=sde$simulate(data=data)
  
  #save plot of the sampled track
  plot_sample=ggplot(sde_sim)+ geom_point(mapping=aes(x, y,col=ID),size= 1,alpha=0.8) +
    scale_color_viridis_d(labels=c("A1","A2","A3","A4","A5","A6"))
  ggsave(paste(paste("check_fe",model_name,"sample",sep="_"),".png",sep=""),plot=plot_sample,width=10,height=5)
  
 
  #fit a new model on the simulated tracks
  sde_check=SDE$new(formulas = formulas,data = sde_sim,type = "sde",response = sde$response(),par0 = par0,map=map)
  sde_check$fit()
  
  #draw coeff from posterior distribution once and for all 
  # to get confidence intervals later
  postcoeff=sde_check$post_coeff(npost)
  all_coeff_re_post=postcoeff$coeff_re
  all_coeff_fe_post=postcoeff$coeff_fe
  
  #get all covariables in the model
  all_vars <- unique(unlist(lapply(formulas, function(formula) all.vars((formula)))))
  
  #Set link function to identity if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(links))) {
      link=(\(x) x)
      links[[var]]=link
    }
  }
  
  #Set xmins and xmaxs to min and max of observed covariates if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(xmins))) {
      xmin=min(data[,var])
      xmins[[var]]=xmin
    }
    if (!(var %in% names(xmaxs))) {
      xmax=max(data[,var])
      xmaxs[[var]]=xmax
    }
  }
  
  #Set xlabels to name of covariates if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(xlabels))) {
      xlabels[[var]]=var
    }
  }
  
  
  #loop over the formulas of each parameter 
  for (j in seq_along(formulas)) {
    
    #formula for the parameter j
    form=formulas[[j]] 
    
    #name of parameter j
    par=pars[j] 
    print(par)
    
    #fixed effects covariates/variables in the formula
    vars=all.vars(form) 
    vars=vars[!(vars %in% c("ID","par0","n"))]
    
    #number of covariates
    ncovs=length(vars) 
      
    #if there is only one fixed effect covariable
    if (ncovs==1){
      
      link=links[[1]]
      xlabel=xlabels[[1]]
      xmax=xmaxs[[vars[1]]]
      xmin=xmins[[vars[1]]]
      
      
      p_check=plot_fe_par_2D(baseline=NULL,sde_check,model_name,par,npoints=200,xmin,xmax,link,xlabel,
                     all_coeff_re_post,all_coeff_fe_post,level,save=FALSE)
      p=plot_fe_par_2D(baseline=NULL,sde,model_name,par,npoints=200,xmin,xmax,link,xlabel,
                       all_coeff_re_post=NULL,all_coeff_fe_post=NULL,level,save=FALSE)
      
      combined_plot <- ggplot() +
        geom_line(data = p_check$data, aes(x =X, y =par, color = "Fixed effects estimate"), size = 1) +
        geom_ribbon(data = p_check$data,aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)+
        geom_line(data = p$data, aes(x =X, y = par, color = "True fixed effects"), size = 1) +
        scale_color_manual(values = c("Fixed effects estimate" = "red", "True fixed effects" = "blue"))+
        xlab(xlabels[[var]])+labs(color = "")
      
      if (!(dir.exists(model_name))) {
        dir.create(model_name)
      }
      
      ggsave(paste(paste("check_fe",model_name,par,var,sep="_"),".png",sep=""),plot=combined_plot,width=10,height=5,path=model_name)
      
      
    }
    
    # if there are two fixed effects covariables
    else if (ncovs==2){
      
      var1=vars[[1]]
      var2=vars[[2]]
      #if the covariates are orthogonal
      if (are_orthogonals(data,var1,var2)) {
        
        list_plots_check=plot_fe_par_2D_ortho(baseline=NULL,sde_check,model_name,par,npoints=200,xmins,xmaxs,links,xlabels,
                             all_coeff_re_post,all_coeff_fe_post,level,save=FALSE)
        list_plots=plot_fe_par_2D_ortho(baseline=NULL,sde,model_name,par,npoints=200,xmins,xmaxs,links,xlabels,
                                        all_coeff_re_post=NULL,all_coeff_fe_post=NULL,level,save=FALSE)
        p1_check=list_plots_check[[1]]
        p2_check=list_plots_check[[2]]
        p1=list_plots[[1]]
        p2=list_plots[[2]]
        
        combined_plot1 <- ggplot() +
          geom_line(data = p1_check$data, aes(x =X, y =par, color = "Fixed effects estimate"), size = 1) +
          geom_ribbon(data = p1_check$data,aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)+
          geom_line(data = p1$data, aes(x =X, y = par, color = "True fixed effects"), size = 1) +
          scale_color_manual(values = c("Fixed effects estimate" = "red", "True fixed effects" = "blue"))+
          xlab(xlabels[[var1]])+labs(color = "")
        
        combined_plot2 <- ggplot() +
          geom_line(data = p2_check$data, aes(x =X, y =par, color = "Fixed effects estimate"), size = 1) +
          geom_ribbon(data = p2_check$data,aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)+
          geom_line(data = p2$data, aes(x =X, y = par, color = "True fixed effects"), size = 1) +
          scale_color_manual(values = c("Fixed effects estimate" = "red", "True fixed effects" = "blue"))+
          xlab(xlabels[[var2]])+labs(color = "")
        
        if (!(dir.exists(model_name))) {
          dir.create(model_name)
        }
        
        ggsave(paste(paste("check_fe",model_name,par,var1,sep="_"),".png",sep=""),
               plot=combined_plot1,width=10,height=5,path=model_name)
        ggsave(paste(paste("check_fe",model_name,par,var2,sep="_"),".png",sep=""),
               plot=combined_plot2,width=10,height=5,path=model_name)
        
        
      }
      
      #else, plot in 3D
      else {
        
        p=plot_fe_par_3D(baseline,sde,model_name,par,npoints=200,xmins=xmins,xmaxs=xmaxs,links=links,
                       xlabels=xlabels,all_coeff_re_post,all_coeff_fe_post,level=0.95,save=FALSE)
      }
    }
  }
}




nearest_shore_point=function(point,coastline) {
  #find the nearest shorepoints on the coastline
  
  # Initialize variables to store the minimum distance and nearest point
  min_distance <- Inf
  nearest_point <- NULL
  
  
  # Iterate over each polygon in coastline_utm
  for (j in 1:length(coastline$geometry)) {
    
    # Get the current polygon
    polygon <- coastline$geometry[[j]]
    
    
    distance <- st_distance(point, polygon)
    
    candidate=st_nearest_points(point,polygon)
    if (!st_is_empty(candidate)) {
      coordinates=st_coordinates(candidate)[2,c("X","Y")]
    
      if (distance < min_distance) {
        min_distance=distance
        nearest_point=coordinates
      }
    }
  }
  return (as.matrix(nearest_point))
}



simulate_constrained_ctcrw=function(z0,ctcrw,data,coastline) {
  
  ###############################################
  
  #z0 : initial point of the simulation
  #ctcrw : fitted ctcrw SDE
  #coastline : coastline data as UTM polygon
  
  #Simulate from a 2D CTCRW whose parameter mu depends
  #on the (normalized) normal vectors between the animal and the shore
  
  ###############################################"
  
  par <- ctcrw$par(new_data = data)
  
  # Loop over dimensions
  n_dim <- length(ctcrw$response())
  if(length(z0) < n_dim) {
    z0 <- rep(z0, n_dim)
  }
  # Initialize vector of simulated observations
  obs <- as.data.frame(matrix(rep(c(NA,NA), nrow(data)),nrow=nrow(data)))
  
  # Loop over IDs
  for(id in seq_along(unique(data$ID))) {
    
    # Get relevant rows of data
    print(id)
    ind <- which(data$ID == unique(data$ID)[id])
    dtimes <- diff(data$time[ind])
    sub_n <- length(ind)
    sub_par <- par[ind,]
    
    # hidden state
    alpha=matrix(c(z0[1],0,z0[2],0),nrow=4,ncol=1,byrow=TRUE)
    
    # Data has two columns for velocity and two for position
    sub_dat=data.frame("x1"=alpha[1],"v1"=alpha[2],"x2"=alpha[3],"v2"=alpha[4])
    
    # Unpack parameters
    tau <- sub_par[, n_dim + 1]
    nu <- sub_par[, n_dim + 2]
    beta <- 1/tau
    sigma <- 2 * nu / sqrt(tau * pi)
    
    # Loop over time steps
    mean <- rep(NA, 2)
    
    for(i in 2:sub_n) {
      
      #last location
      z=as.numeric(sub_dat[i-1,c("x1","x2")])
      #compute normal vector
      nearest_point=nearest_shore_point(st_point(1000*z),coastline)/1000
      nx=z[1]-nearest_point[1]
      ny=z[2]-nearest_point[2]
      Dshore=sqrt(nx^2+ny^2)
      normalized_nx=nx/Dshore
      normalized_ny=ny/Dshore
      
      #adjust covariate data
      new_data=data[i-1,]
      new_data$normalized_nx=normalized_nx
      new_data$normalized_ny=normalized_ny
      
      #get value of mu
      mu=ctcrw$par(new_data=new_data)[1,c("mu1","mu2")]
      mu=matrix(mu,ncol=1,nrow=2,byrow=TRUE)
      
      
      #covariance matrix
      var_xi=sigma[i-1]^2/beta[i-1]^2*(dtimes[i-1]+2/beta[i-1]*(exp(-beta[i-1]*dtimes[i-1])-1)+
                                         1/2/beta[i-1]*(1-exp(-2*beta[i-1]*dtimes[i-1])))
      var_zeta=sigma[i-1]^2/2/beta[i-1]*(1-exp(-2*beta[i-1]*dtimes[i-1]))
      cov=sigma[i-1]^2/2/beta[i-1]^2*(1+exp(-2*beta[i-1]*dtimes[i-1])-2*exp(-beta[i-1]*dtimes[i-1]))
      Qi=matrix(c(var_xi,cov,cov,var_zeta),nrow=2,byrow=TRUE)
      Qi_tilde=bdiag(replicate(2,Qi,simplify=FALSE))
      
      #drift part
      Bi=matrix(c(dtimes[i-1],0),nrow=2,ncol=1,byrow=TRUE)
      Bi_tilde=bdiag(replicate(2,Bi,simplify=FALSE))
      
      #random part
      eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi_tilde)
      
      #hidden states link matrix
      Ti=matrix(c(1,1/beta[i-1]*(1-exp(-beta[i-1]*dtimes[i-1])),0,exp(-beta[i-1]*dtimes[i-1])),nrow=2,byrow=TRUE)
      Ti_tilde=bdiag(replicate(2,Ti,simplify=FALSE))
      
      #state obs link matrix
      Z=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE)
      Z_tilde=bdiag(replicate(2,Z,simplify=FALSE))
      
      
      alpha=Ti_tilde%*%alpha+Bi_tilde%*%mu+eta
      sub_dat=rbind(sub_dat,c(alpha[1,1],alpha[2,1]+mu[1],alpha[3,1],alpha[4,1]+mu[2]))
      
    }
    # Only return location (and not velocity)
    sub_obs <- sub_dat[,c("x1","x2")]
    # Update observation vector
    obs[ind,] <- sub_obs
  }

  # Add simulated variable to data frame
  data[,ctcrw$response()] <- obs
  return (data)
}


#' @description Compute signed angle in [-pi,pi] that rotates first vector into second vector
#' 
#' @param point the location we want to tell if it is on land or not
#' @param land sf object (list of polygons) defining the land
#' @return boolean TRUE if it is on land, FALSE otherwise
#' 
signed_angle <- function(u, v) {
  # u, v: matrices of bivariate vectors with 2 columns
  u <- matrix(u, ncol = 2)
  v <- matrix(v, ncol = 2)
  if (nrow(u) != nrow(v)) stop("u and v must have the same number of 
                                  rows")
  result <- as.numeric(atan2(v[,2], v[,1]) - atan2(u[,2], u[,1]))
  ind1 <- which(result > pi)
  ind2 <- which(result <= -pi)
  result[ind1] <- result[ind1] - 2*pi
  result[ind2] <- result[ind2] + 2*pi
  return(result) 
} 
