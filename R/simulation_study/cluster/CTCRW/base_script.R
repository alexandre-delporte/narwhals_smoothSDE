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
# Script Description: Base script for simulation study with CTCRW.. 
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
par_dir=here("R","simulation_study","cluster","CTCRW",domain_name)
set_up_file=paste("set_up_",domain_name,".R",sep="")
source(file.path(par_dir,set_up_file))     #get set up for simulation study
source(file.path(here("R","simulation_study","CVM_functions.R")))  #get functions to simulate trajectories


seed= 1
set.seed(seed)




# Generate samples -----------------------------

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


data=foreach (i=1:N_ID_HIGH,.combine='rbind',.packages=c("progress","MASS","sf","mgcv","smoothSDE")) %dopar% {
  
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
if (!(dir.exists(file.path(par_dir,paste("results_",hyper_params_file_name,sep=""))))) {
  dir.create(file.path(par_dir,paste("results_",hyper_params_file_name,sep="")))
}

ggsave(path=file.path(par_dir,paste("results_",hyper_params_file_name,sep="")),
       filename=paste("simulated_trajectories_","seed",seed,"_",
                      hyper_params_file_name,".png",sep=""),plot=plot,width=10,height=5)

if (count>0) {
  stop("Stop : at least one trajectory reached the shore.")
}


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
clusterExport(cl, varlist = c("formulas_ctcrw", "new_lambda_ctcrw", 
                              "PAR0", "SIGMA_OBS_LOW","SIGMA_OBS_HIGH","N_ID_LOW"))

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
      par0 = c(0,0,PAR0),
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


write_estimates_csv=function(results,model_type="ctcrw") {
  
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
      
      # add mean obtained from samples with joint covariance matrix
      post_coeff=try(model$post_coeff(n_post=10000))
      if (inherits(post_coeff, "try-error")) {
        message(paste("Invalid joint covariance matrix for ", model_name,'\n', model))
        next
      }
      
      post_par=list(
        "tau"=exp(post_coeff$coeff_fe[,"tau.(Intercept)"]),
        "nu"=exp(post_coeff$coeff_fe[,"nu.(Intercept)"]),
        "sigma_tau"=1/sqrt(exp(post_coeff$log_lambda[,"tau.s(ID)"])),
        "sigma_nu"=1/sqrt(exp(post_coeff$log_lambda[,"nu.s(ID)"])))
      
      mean_par<-lapply(post_par,mean)
      sd_par<-lapply(post_par,sd)
      post_par_df<- data.frame("coeff_name"=names(post_par),"estimate"=unlist(mean_par),
                               "std"=unlist(sd_par))
      rownames(post_par_df) <- NULL
      
      coeffs_df=rbind(coeffs_df,post_par_df)
      
      
      # path for csv file
      output_file <- file.path(par_dir, paste0("results_",hyper_params_file_name),
                               paste0("estimates_", model_name, "_seed", seed, ".csv"))
      
      # Wwite the csv file
      write.csv(coeffs_df, file = output_file, row.names = FALSE)
    }
  }
  
}

write_estimates_csv(ctcrw_results,"ctcrw")
