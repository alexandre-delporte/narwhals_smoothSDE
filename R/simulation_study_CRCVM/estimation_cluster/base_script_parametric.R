# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2025 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2025-01-15
#
# Script Name:    
#
# Script Description:
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(here)               # To get root of git repo
library(smoothSDE)          #to compute nearest shore point
library(foreach)            #foreach loop
library(doParallel)         #parallel computing



domain_name="fjords"
type="parametric"
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
  data_shore=res$shore
  data_sim=cbind(data_sim,data_shore)
  data_sim
}

#stop cluster
stopCluster(cl)

#rename Dshore
colnames(data)[7]="DistanceShore"

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
if (!(dir.exists(file.path(par_dir,paste("results_","parametric_",hyper_params_file_name,sep=""))))) {
  dir.create(file.path(par_dir,paste("results_","parametric_",hyper_params_file_name,sep="")))
}

ggsave(path=file.path(par_dir,paste("results_","parametric_",hyper_params_file_name,sep="")),
       filename=paste("simulated_trajectories_","seed",seed,"_",
                      hyper_params_file_name,".png",sep=""),plot=plot,width=10,height=5)

if (count>0) {
  stop("Stop : at least one trajectory reached the shore.")
}


all_data=list("data_lf_he"=data_lf_he,"data_hf_le"=data_hf_le,"data_lf_le"=data_lf_le,
              "data_hf_he"=data_hf_he)


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
        message(paste("Skipping invalid model ", model_name,'\n', model))
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
      
      # path for csv file
      output_file <- file.path(par_dir, paste0("results_","parametric_",hyper_params_file_name),
                               paste0("estimates_", model_name, "_seed", seed, ".csv"))
      
      # Wwite the csv file
      write.csv(coeffs_df, file = output_file, row.names = FALSE)
    }
  }
  
}


estimate_crcvm_parallel=function(data_list) {
  # Set up the parallel backend
  n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to the cluster
  clusterExport(cl, varlist = c("formulas_crcvm", "new_lambda_crcvm","PAR0", "A","B","D0","D1",
                                "SIGMA_D","SIGMA_THETA","SIGMA_OBS_LOW",
                                "SIGMA_OBS_HIGH", "N_ID_LOW"))
  
  # Parallel execution
  crcvm_results <- foreach(
    data_name = names(data_list),
    .combine = rbind,
    .packages = c("smoothSDE")
  ) %dopar% {
    dataset <- data_list[[data_name]]
    
    error_type=strsplit(data_name,"_")[[1]][3]
    sigma_obs=ifelse(error_type=="le",SIGMA_OBS_LOW,SIGMA_OBS_HIGH)
    H_high=array(rep(sigma_obs^2*diag(2),length(dataset$time)),dim=c(2,2,length(dataset$time)))
    
    # High number of IDs
    fit_crcvm<- function(data,H) {
      sde <- SDE$new(formulas = formulas_crcvm,data = data,type = "CRCVM_SSM",
                     response = c("y1","y2"),
                     par0 =c(PAR0,A,B,D0,D1,SIGMA_D,SIGMA_THETA),
                     fixpar=c("a","b","D0","D1","sigma_D","sigma_theta"),
                     other_data=list("H"=H))
      
      sde$update_lambda(new_lambda_crcvm)
      sde$fit(method = "BFGS")
      
      return (sde)
    }
    crcvm_high<-try(fit_crcvm(dataset,H_high))
    if (!(inherits(crcvm_high, "SDE"))) {
      message(paste("Error in fit with high nb of ID for dataset:", data_name, "\nError message:", crcvm_high))
    }
    
    
    # Low number of IDs
    dataset_low <- dataset[dataset$ID %in% 1:N_ID_LOW, ]
    H_low=array(rep(sigma_obs^2*diag(2),length(dataset_low$time)),dim=c(2,2,length(dataset_low$time)))
    
    crcvm_low<-try(fit_crcvm(dataset_low,H_low))
    
    if (!(inherits(crcvm_low, "SDE"))) {
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
  
  return(crcvm_results)
}


# Estimate from simulated data with true DistanceShore and theta  --------------------------------------------

cat("Estimating CRCVM parameters with true DistanceShore and theta...\n")

# Parameters and formulas
formulas_crcvm <- list(tau=~s(ID,bs="re"),nu=~s(ID,bs="re"),a=~1,b=~1,
                       D0=~1,D1=~1,sigma_D=~1,sigma_theta=~1)
new_lambda_crcvm=c(1/SIGMA_TAU_0^2,1/SIGMA_NU_0^2)

crcvm_true_results<-estimate_crcvm_parallel(all_data)

cat("Saving estimates in csv file...\n")

# Write estimates parameters in csv file
write_estimates_csv(crcvm_true_results,"crcvm_true")


# Estimate from simulated data with observed DistanceShore and theta

# Function to calculate rolling average
rolling_average <- function(data, window) {
  n <- nrow(data)
  result <- matrix(NA, n, ncol(data))  # Initialize result matrix with NA
  colnames(result) <- colnames(data)
  
  for (i in seq_len(n)) {
    # Define the range for the current window
    start <- max(1, i - floor(window / 2))
    end <- min(n, i + floor(window / 2))
    result[i, ] <- colMeans(data[start:end, , drop = FALSE])  # Compute column means
  }
  return(result)
}


add_covs_parallel <- function(data,add="all",window=1,n_cores = parallel::detectCores() - 1) {
  # Setup parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Split data by ID
  ids <- unique(data$ID)
  
  # Process each ID in parallel
  results <- foreach(id = ids, .combine = rbind, .packages = c("sf"),
                     .export = c("nearest_boundary_points", "is_in_land", "signed_angle","border","rolling_average")  ) %dopar% {
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
                                   
                                   
                                   vexp_avg <- rolling_average(vexp_df, window)
                                   
                                   # Nearest points on shore
                                   sub_data <- cbind(sub_data, nearest_boundary_points(as.matrix(sub_data[, c("y1", "y2")]), border))
                                   
                                   # Normal vectors
                                   normal <- as.matrix(sub_data[2:n_sub, c("y1", "y2")] - sub_data[2:n_sub, c("p1", "p2")])
                                   
                                   # Angle between velocity and normal vector
                                   theta_coast <- signed_angle(normal, vexp_avg)
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
                                   
                                   Dshore_avg=as.vector(rolling_average(matrix(Dshore,ncol=1),window))
                                   
                                   # Calculate Ishore
                                   Ishore_avg <- 1/Dshore_avg
                                   
                                   # Add columns to sub_data
                                   if (add=="all") {
                                     sub_data$theta <- theta_coast
                                     sub_data$DistanceShore <- Dshore_avg
                                     sub_data$Ishore <- Ishore_avg
                                   }
                                   else if (add=="theta") {
                                     sub_data$theta <- theta_coast
                                   }
                                   else if (add=="DistanceShore") {
                                     sub_data$DistanceShore <- Dshore
                                     sub_data$Ishore <- Ishore
                                   }
                                   else {
                                     stop("Options for argument add are \"all\",\"DistanceShore\" and \"theta\"")
                                   }
                                   
                                   return(sub_data)
                                 }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(results)
}
#replace true DistanceShore and theta by rolling average of observed DistanceShore and theta 
data_lf_he_obs=add_covs_parallel(data_lf_he)
data_hf_le_obs=add_covs_parallel(data_hf_le)
data_lf_le_obs=add_covs_parallel(data_lf_le)
data_hf_he_obs=add_covs_parallel(data_hf_he)

all_data_obs=list("data_lf_he"=data_lf_he_obs,"data_hf_le"=data_hf_le_obs,"data_lf_le"=data_lf_le_obs,
              "data_hf_he"=data_hf_he_obs)


cat("Estimating CRCVM parameters with observed DistanceShore and theta...\n")

crcvm_obs_results<-estimate_crcvm_parallel(all_data_obs)


cat("Saving estimates in csv file...\n")

# Write estimates parameters in csv file
write_estimates_csv(crcvm_obs_results,"crcvm_obs")


#replace true DistanceShore and theta by observed DistanceShore and theta 
data_lf_he_avg=add_covs_parallel(data_lf_he,window=5)
data_hf_le_avg=add_covs_parallel(data_hf_le,window=5)
data_lf_le_avg=add_covs_parallel(data_lf_le,window=5)
data_hf_he_avg=add_covs_parallel(data_hf_he,window=5)


all_data_avg=list("data_lf_he"=data_lf_he_avg,"data_hf_le"=data_hf_le_avg,"data_lf_le"=data_lf_le_avg,
                  "data_hf_he"=data_hf_he_avg)


cat("Estimating CRCVM parameters with rolling average of observed DistanceShore and theta...\n")

crcvm_avg_results<-estimate_crcvm_parallel(all_data_avg)


cat("Saving estimates in csv file...\n")

# Write estimates parameters in csv file
write_estimates_csv(crcvm_avg_results,"crcvm_avg")

