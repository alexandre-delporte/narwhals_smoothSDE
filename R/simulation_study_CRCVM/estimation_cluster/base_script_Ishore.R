# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-03
#
# Script Name:    base_script_Ishore.R
#
# Script Description: Base script for simulation study with Ishore covariate.. 
# Scripts to execute on the clusters for parameter estimation are generated from this script.
#
#
# SETUP --------------------------------------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(here)               # To get root of git repo
library(smoothSDE)          #to compute nearest shore point
library(foreach)            #foreach loop
library(doParallel)         #parallel computing



domain_name="fjords"
type="Ishore"
par_dir=here("R","simulation_study_CRCVM","estimation_cluster",domain_name)
set_up_file=paste("set_up_",domain_name,"_",type,".R",sep="")
source(file.path(par_dir,set_up_file))     #get set up for simulation study
source(file.path(here("R","simulation_study_CRCVM","CVM_functions.R")))  #get functions to simulate trajectories


seed= 1
set.seed(seed)




# Generate samples -----------------------------

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


data=foreach (i=1:N_ID_HIGH,.combine='rbind',.packages=c("progress","MASS","sf","mgcv")) %dopar% {
  
  set.seed((seed-1)*N_ID_HIGH+i)
  
  #constant nu
  fnu_constant=function(cov_data) {
    return (exp(true_log_nu[i]))
  }
  
  ftau_constant=function(cov_data) {
    return (exp(true_log_tau[i]))
  }
  res=sim_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,
                log_sigma_obs=NULL,v0=v0,x0=x0[i,],times=times,land=border,verbose=FALSE)
  
  data_sim=res$sim
  data_sim$ID=factor(rep(i,length(data_sim$y1)))
  data_sim
}

#stop cluster
stopCluster(cl)



# Points that reached land 
count=0
for (id in unique(data$ID)) {
  sub_data=data[data$ID==id,]
  if (nrow(sub_data) < n_obs) {
    count=count+1
    cat("ID",id,"reached land","\n",sep=" ")
  }
}

cat(count/N_ID_HIGH*100,"percent of the samples reached land \n")

if (count>0) {
  stop("At least one trajectory reached land.")
}


# add noise
low_noise=rmvn(n_obs,rep(0,2),diag(rep(SIGMA_OBS_LOW^2,2)))
high_noise=rmvn(n_obs,rep(0,2),diag(rep(SIGMA_OBS_HIGH^2,2)))
observed_data_le=data
observed_data_le[,c("y1","y2")]=data[,c("y1","y2")]+low_noise
observed_data_he=data
observed_data_he[,c("y1","y2")]=data[,c("y1","y2")]+high_noise

# subsample
data_lf_he=observed_data_he[seq(1,length(data$time),by=BY_LF),]
data_lf_le=observed_data_le[seq(1,length(data$time),by=BY_LF),]
data_hf_le=observed_data_le[seq(1,length(data$time),by=BY_HF),]
data_hf_he=observed_data_he[seq(1,length(data$time),by=BY_HF),]


# Save plot of the trajectories

plot=ggplot()+geom_sf(data=border$geometry,fill="lightgrey")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_point(data=data,mapping=aes(y1,y2,col=ID),size=0.1)+
  geom_path(data=data,mapping=aes(y1,y2,color=ID),size=0.1)+
  geom_point(data = data%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

hyper_params_file_name=sapply(strsplit(hyperparams_file,"\\."),'[',1)

#create directory to save the results
if (!(dir.exists(file.path(par_dir,paste("results_","Ishore_",hyper_params_file_name,sep=""))))) {
  dir.create(file.path(par_dir,paste("results_","Ishore_",hyper_params_file_name,sep="")))
}

ggsave(path=file.path(par_dir,paste("results_","Ishore_",hyper_params_file_name,sep="")),
       filename=paste("simulated_trajectories_","seed",seed,"_",
                      hyper_params_file_name,".png",sep=""),plot=plot,width=10,height=5)

if (count>0) {
  stop("Stop : at least one trajectory reached the shore.")
}

# Add covariates to simulation -------------------------------------------------

signed_angle <- function(u, v) {
  #Compute signed angle in [-pi,pi] that rotates first vector into second vector 
  # as in 
  # https://math.stackexchange.com/questions/529555/signed-angle-between-2-vectors
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


add_covs_parallel <- function(data, n_cores = parallel::detectCores() - 1) {
  # Setup parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Split data by ID
  ids <- unique(data$ID)
  
  # Process each ID in parallel
  results <- foreach(id = ids, .combine = rbind, .packages = c("sf"),
                     .export = c("nearest_boundary_points", "is_in_land", "signed_angle","border",
                                 "D_LOW","D_UP")  ) %dopar% {
                                   # Filter data for the current ID
                                   sub_data <- data[data$ID == id, ]
                                   n_sub <- nrow(sub_data)
                                   
                                   # Time steps
                                   dtimes <- sub_data[2:n_sub, "time"] - sub_data[1:(n_sub - 1), "time"]
                                   
                                   # Step lengths
                                   dx <- sub_data[2:n_sub, "y1"] - sub_data[1:(n_sub - 1), "y1"]
                                   dy <- sub_data[2:n_sub, "y2"] - sub_data[1:(n_sub - 1), "y2"]
                                   
                                   # Empirical velocity
                                   vexp_df <- cbind(dx / dtimes, dy / dtimes)
                                   
                                   # Nearest points on shore
                                   sub_data <- cbind(sub_data, nearest_boundary_points(as.matrix(sub_data[, c("y1", "y2")]), border))
                                   
                                   # Normal vectors
                                   normal <- as.matrix(sub_data[2:n_sub, c("y1", "y2")] - sub_data[2:n_sub, c("p1", "p2")])
                                   
                                   # Angle between velocity and normal vector
                                   theta_coast <- signed_angle(normal, vexp_df)
                                   theta_coast <- c(theta_coast, 1)  # Adjust length
                                   
                                   # Initialize DistanceShore
                                   Dshore <- rep(0, n_sub)
                                   
                                   for (i in 1:n_sub) {
                                     y1 <- sub_data[i, "y1"]
                                     y2 <- sub_data[i, "y2"]
                                     if (!is_in_land(st_point(c(y1, y2)), border)) {
                                       p1 <- sub_data[i, "p1"]
                                       p2 <- sub_data[i, "p2"]
                                       Dshore[i] <- sqrt((y1 - p1)^2 + (y2 - p2)^2)
                                     }
                                   }
                                   
                                   # Calculate Ishore
                                   Ishore <- ifelse(Dshore < D_LOW, 1 / D_LOW, ifelse(Dshore > D_UP, 0, 1 / Dshore))
                                   
                                   # Add columns to sub_data
                                   sub_data$theta <- theta_coast
                                   sub_data$DistanceShore <- Dshore
                                   sub_data$Ishore <- Ishore
                                   
                                   return(sub_data)
                                 }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(results)
}

data_lf_he=add_covs_parallel(data_lf_he)
data_hf_le=add_covs_parallel(data_hf_le)
data_lf_le=add_covs_parallel(data_lf_le)
data_hf_he=add_covs_parallel(data_hf_he)

all_data=list("data_lf_he"=data_lf_he,"data_hf_le"=data_hf_le,"data_lf_le"=data_lf_le,
              "data_hf_he"=data_hf_he)


# Estimate from simulated data with CTCRW -----------------------------------------------------------

cat("\n Estimating CTCRW parameters...\n")


# Set up the parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Parameters and formulas
formulas_ctcrw <- list(mu1 = ~1, mu2 = ~1, tau = ~s(ID, bs = "re"), nu = ~s(ID, bs = "re"))
new_lambda_ctcrw <- c(1 / SIGMA_TAU_0^2, 1 / SIGMA_NU_0^2)

# Export necessary objects to the cluster
clusterExport(cl, varlist = c("formulas_ctcrw", "new_lambda_ctcrw", "PAR0", "SIGMA_OBS_LOW","SIGMA_OBS_HIGH","N_ID_LOW"))

# Parallel execution
ctcrw_results <- foreach(
  data_name = names(all_data),
  .combine = rbind,
  .packages = c("smoothSDE")
) %dopar% {
  dataset <- all_data[[data_name]]
  
  error_type=strsplit(data_name,"_")[[1]][3]
  sigma_obs=ifelse(error_type=="le",SIGMA_OBS_LOW,SIGMA_OBS_HIGH)
  H_high=array(rep(sigma_obs^2*diag(2),length(dataset$time)),dim=c(2,2,length(dataset$time)))
  
  # High number of IDs
  fit_ctcrw<- function(data,H) {
    
    sde<- SDE$new(
      formulas = formulas_ctcrw,
      data = data,
      type = "CTCRW",
      response = c("y1", "y2"),
      par0 = PAR0[1:4],
      fixpar = c("mu1", "mu2"),
      other_data = list("H" = H)
    )
    sde$update_lambda(new_lambda_ctcrw)
    sde$fit(method = "BFGS")
    
    return(sde)
  }
  ctcrw_high<-try(fit_ctcrw(dataset,H_high),silent=TRUE)
  
  if (inherits(ctcrw_high, "try-error")) {
    message(paste("Error in fit with high nb of ID for dataset:", data_name, "\nError message:", ctcrw_high))
  }
  
  
  # Low number of IDs
  dataset_low <- dataset[dataset$ID %in% 1:N_ID_LOW, ]
  H_low=array(rep(sigma_obs^2*diag(2),length(dataset_low$time)),dim=c(2,2,length(dataset_low$time)))
  
  ctcrw_low<-try(fit_ctcrw(dataset_low,H_low))
  
  if (inherits(ctcrw_low, "try-error")) {
    message(paste("Error in fit with low nb of ID for dataset:", data_name, "\nError message:", ctcrw_low))
  }
  
  # Combine results
  list(
    data_name = data_name,
    fit_high = ctcrw_high,
    fit_low = ctcrw_low
  )
}

stopCluster(cl)

cat("Saving estimates in csv file...\n")

# Write estimated parameters in csv files -----------------------------


estimate_parametric_omega=function(model,fixed_a=1,probs=c(0.7,0.95)) {
  
  data=model$data()
  
  Dshore=seq(from=quantile(data$DistanceShore,probs[1]),
             to=quantile(data$DistanceShore,probs[2]),length.out=30)
  Ishore=1/Dshore
  theta=seq(from=-pi,to=pi,length.out=30)
  grid<- as.data.frame(expand.grid(Ishore,theta))
  colnames(grid) <- c("Ishore","theta")
  grid$ID=unique(data$ID)[1]
  grid$Dshore=1/grid$Ishore
  
  #get model matrices
  mats=model$make_mat(new_data=grid)
  X_fe=mats$X_fe
  X_re=mats$X_re
  
  par_mat=model$par(new_data=grid,X_fe=X_fe,X_re=X_re)
  
  #matrix of values for surface plot
  z=matrix(par_mat[,"omega"],30,30)
  
  data=grid
  data$z=as.vector(z)
  
  
  
  data_copy=data
  colnames(data_copy)=c("Ishore","theta","ID","DistanceShore","z")
  
  fit_optim <- optim(par = log(c(0.5,0.1,2,0.5,pi/6)),
                     fn = function(params) {
                       
                       # Compute predicted z values
                       z_pred <- fomega_cubic(data_copy,fixed_a,exp(params[1]),
                                              exp(params[2]),exp(params[3]),exp(params[4]),exp(params[5]))
                       
                       # Return sum of squared residuals
                       sum((data_copy$z - z_pred)^2)
                     })
  
  
  est_par=c(fixed_a,exp(fit_optim$par))
  names(est_par)=c("a","b","D0","D1","sigma_D","sigma_theta")
  
  return (est_par)
  
}

write_estimates_csv=function(results,model_type) {
  
  n=dim(results)[1]
  p=dim(results)[2]
  
  for (i in 1:n) {
    
    data_name=results[i,1][[1]]
    settings=strsplit(data_name,"_")[[1]][2:3]
    
    for (j in 1:(p-1)) {
      
      chain=ifelse(j==1,"hID","lID")
      model=results[i,j+1][[1]]
    
      model_name=paste(model_type,settings[1],settings[2],chain,sep="_")
      
      # Check if model is a valid SDE and skip if not
      if (!inherits(model, "SDE") || inherits(model, "try-error")) {
        message(paste("Skipping invalid model ", model_name))
        next
      }
      
      #re and fe coeffs
      coeffs=rbind(model$coeff_re(),model$coeff_fe()) 
      
      coeff_names=rownames(coeffs)
      coeff_values=as.numeric(coeffs)
      
      #standard deviation of random effects as smoothing penalties
      log_lambda=log(model$lambda())
      log_lambda_names=rownames(log_lambda)
      log_lambda_values=as.numeric(log_lambda)
      
      #measurement error
      log_sigma_obs_value=as.list(model$tmb_rep(),what="Est")$log_sigma_obs
      
      #standard errors
      all_std=as.list(model$tmb_rep(),what="Std")
      std_coeffs=rbind(all_std$coeff_re,all_std$coeff_fe)
      std_log_lambda=all_std$log_lambda
      std_log_sigma_obs=all_std$log_sigma_ob
      
      #create dataframe
      coeffs_df=data.frame("coeff_name"=factor(c(coeff_names,log_lambda_names,"log_sigma_obs")),
                           "estimate"=c(coeff_values,log_lambda_values,log_sigma_obs_value),
                           "std"=c(std_coeffs,std_log_lambda,std_log_sigma_obs))
      
      if (model_type=="crcvm") {
        
        est_par=estimate_parametric_omega(model,fixed_a=A)
        
        new_coeffs_df=data.frame("coeff_name"=factor(names(est_par)),"estimate"=as.numeric(est_par),
                                 "std"=rep(NA,length(est_par)))
        
        coeffs_df=rbind(coeffs_df,new_coeffs_df)
      }
      
      # path for csv file
      output_file <- file.path(par_dir, paste0("results_","Ishore_",hyper_params_file_name),
                               paste0("estimates_", model_name, "_seed", seed, ".csv"))
      
      # Wwite the csv file
      write.csv(coeffs_df, file = output_file, row.names = FALSE)
    }
  }
  
}

write_estimates_csv(ctcrw_results,"ctcrw")



# Estimate from simulated data with CRCVM  --------------------------------------------

cat("Estimating CRCVM parameters...\n")

# Set up the parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Parameters and formulas
formulas_crcvm <- list(mu1=~1,mu2=~1,tau=~s(ID,bs="re"),
                       nu=~s(ID,bs="re"),omega=~te(theta,Ishore,k=SP_DF,bs="cs"))
new_lambda_crcvm=c(1/SIGMA_TAU_0^2,1/SIGMA_NU_0^2,lambda_splines)
new_map_crcvm=list("log_lambda"=factor(c(1,2,NA,NA)))

# Export necessary objects to the cluster
clusterExport(cl, varlist = c("formulas_crcvm", "new_lambda_crcvm", "new_map_crcvm", "PAR0", "SIGMA_OBS_LOW",
                              "SIGMA_OBS_HIGH", "N_ID_LOW", "SP_DF"))

# Parallel execution
crcvm_results <- foreach(
  data_name = names(all_data),
  .combine = rbind,
  .packages = c("smoothSDE")
) %dopar% {
  dataset <- all_data[[data_name]]
  
  error_type=strsplit(data_name,"_")[[1]][3]
  sigma_obs=ifelse(error_type=="le",SIGMA_OBS_LOW,SIGMA_OBS_HIGH)
  H_high=array(rep(sigma_obs^2*diag(2),length(dataset$time)),dim=c(2,2,length(dataset$time)))
  
  # High number of IDs
  fit_crcvm<- function(data,H) {
    sde <- SDE$new(formulas = formulas_crcvm,data = data,type = "RACVM_SSM",
                   response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                   other_data=list("H"=H))
    
    sde$update_lambda(new_lambda_crcvm)
    sde$update_map(new_map_crcvm)
    sde$fit(method = "BFGS")
    
    return (sde)
  }
  crcvm_high<-try(fit_crcvm(dataset,H_high))
  if (inherits(crcvm_high, "try-error")) {
    message(paste("Error in fit with high nb of ID for dataset:", data_name, "\nError message:", crcvm_high))
  }
  
  
  # Low number of IDs
  dataset_low <- dataset[dataset$ID %in% 1:N_ID_LOW, ]
  H_low=array(rep(sigma_obs^2*diag(2),length(dataset_low$time)),dim=c(2,2,length(dataset_low$time)))
  
  crcvm_low<-try(fit_crcvm(dataset_low,H_low))
  
  if (inherits(crcvm_low, "try-error")) {
    message(paste("Error in fit with low nb of ID for dataset:", data_name, "\nError message:", crcvm_low))
  }
  
  # Combine results
  list(
    data_name = data_name,
    fit_high = crcvm_high,
    fit_low = crcvm_low
  )
}

stopCluster(cl)

cat("Saving estimates in csv file...\n")

# Write estimates parameters in csv file
write_estimates_csv(crcvm_results,"crcvm")

# Save spline surface estimates

write_surface=function(results,model_type) {
  
  n=dim(results)[1]
  p=dim(results)[2]
  
  for (i in 1:n) {
    
    data_name=results[i,1][[1]]
    settings=strsplit(data_name,"_")[[1]][2:3]
    
    for (j in 1:(p-1)) {
      
      chain=ifelse(j==1,"hID","lID")
      model=results[i,j+1][[1]]
      model_name=paste(model_type,settings[1],settings[2],chain,sep="_")
      
      # Check if model is a valid SDE and skip if not
      if (!inherits(model, "SDE") || inherits(model, "try-error")) {
        message(paste("Skipping invalid model ", model_name))
        next
      }
      
      xmin=list("Ishore"=quantile(model$data()$Ishore,0.25))
      xmax=list("Ishore"=quantile(model$data()$Ishore,0.95))
      plots=model$get_all_plots(link=list("Ishore"=(\(x) 1/x)),xmin=xmin,xmax=xmax)
      
      surface=plots$fe_omega_theta_Ishore
      
      Dshore_values = as.numeric(surface$x$data[[1]]$x)
      theta_values = as.numeric(surface$x$data[[1]]$y)
      omega_values = as.numeric(surface$x$data[[1]]$z)
      
      data_surface=expand.grid(Dshore = Dshore_values, theta = theta_values)
      data_surface$omega=as.vector(omega_values)
      
      # path for csv file
      output_file <- file.path(par_dir, paste0("results_","Ishore_", hyper_params_file_name),
                               paste0("surface_", model_name, "_seed", seed, ".csv"))
      
      
      write.csv(data_surface,output_file,
                row.names=FALSE)
    }
  }
}

write_surface(crcvm_results,"crcvm")

