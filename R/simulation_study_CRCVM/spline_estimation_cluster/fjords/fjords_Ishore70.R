# SETUP --------------------------------------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(here)               # To get root of git repo
library(smoothSDE)          #to compute nearest shore point
library(foreach)            #foreach loop
library(doParallel)         #parallel computing


domain_name="fjords"
par_dir=here("R","simulation_study_CRCVM","spline_estimation_cluster",domain_name)
set_up_file=paste("set_up_",domain_name,".R",sep="")
source(file.path(par_dir,set_up_file))     #get set up for simulation study
source(file.path(here("R","simulation_study_CRCVM","CVM_functions.R")))  #get functions to simulate trajectories


seed= 70
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

cat(count/N_ID_HIGH*100,"percent of the samples reached land")

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

plot=ggplot()+geom_sf(data=border$geometry,fill="lightgrey",border="black")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_point(data=data,mapping=aes(y1,y2,col=ID),size=0.1)+
  geom_path(data=data,mapping=aes(y1,y2,color=ID),size=0.1)+
  geom_point(data = data%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

ggsave(filename=paste("simulated_trajectories_","seed",seed,
                      sapply(strsplit(hyperparams_file,"\\."),'[',1),".png",sep=""),plot=plot,width=10,height=5)

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


add_covs=function(data) {
  
  #Returns new data frame with three more columns "theta", "DistanceShore" and "Ishore"
  
  new_data=data
  #observed deviation angles
  n=length(new_data$time)
  
  new_data$theta=rep(NA,n)
  
  for (id in unique(new_data$ID)) {
    
    #filter data with specific id
    sub_ind=(new_data$ID==id)
    sub_data=new_data[sub_ind,]
    n_sub=length(sub_data$time)
    
    #time steps
    dtimes=sub_data[2:n_sub,"time"]-sub_data[1:(n_sub-1),"time"]
    
    #step lengths
    dx=sub_data[2:n_sub,"y1"]-sub_data[1:(n_sub-1),"y1"]
    dy=sub_data[2:n_sub,"y2"]-sub_data[1:(n_sub-1),"y2"]
    
    #matrix of empirical velocity
    vexp_df=cbind(dx/dtimes,dy/dtimes)
    
    #nearest points on shore
    sub_data=cbind(sub_data,nearest_boundary_points(as.matrix(sub_data[,c("y1","y2")]),border))
    
    #matrix of normal vectors
    normal=as.matrix(sub_data[2:n_sub,c("y1","y2")]-sub_data[2:n_sub,c("p1","p2")])
    
    #angle between velocity and normal vector
    theta_coast=signed_angle(normal,vexp_df)
    
    #adjust lengths 
    theta_coast=c(theta_coast,1) 
    
    new_data[sub_ind,"theta"]=theta_coast
    
    Dshore=rep(0,n_sub)
    for (i in 1:n_sub) {
      y1=sub_data[i,"y1"]
      y2=sub_data[i,"y2"]
      if (!(is_in_land(st_point(c(y1,y2)),border))) {
        p1=sub_data[i,"p1"]
        p2=sub_data[i,"p2"]
        Dshore[i]=sqrt((y1-p1)^2+(y2-p2)^2)
      }
    }
    new_data[sub_ind,"DistanceShore"]=Dshore
    new_data[sub_ind,"Ishore"]=ifelse(Dshore<D_LOW,1/D_LOW,ifelse(Dshore>D_UP,0,1/Dshore))
  }
  
  return(new_data)
}

data_lf_he=add_covs(data_lf_he)
data_hf_le=add_covs(data_hf_le)
data_lf_le=add_covs(data_lf_le)
data_hf_he=add_covs(data_hf_he)



# Estimate from simulated data with CTCRW -----------------------------------------------------------


## Low frequency high error -------------

formulas <- list(mu1=~1,mu2=~1,tau=~s(ID,bs="re"),
                 nu=~s(ID,bs="re"))


### High number of ID
ctcrw_lf_he_hID<- SDE$new(formulas = formulas,data = data_lf_he,type = "CTCRW",
                response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
new_lambda=c(1/SIGMA_TAU_0^2,1/SIGMA_NU_0^2)
ctcrw_lf_he_hID$update_lambda(new_lambda)

#fit
ctcrw_lf_he_hID$fit(method="BFGS")

### Low number of ID
ctcrw_lf_he_lID<- SDE$new(formulas = formulas,data = data_lf_he[data_lf_he$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                             response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                             other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
ctcrw_lf_he_lID$update_lambda(new_lambda)

#fit
ctcrw_lf_he_lID$fit(method="BFGS")

## Low frequency low error ------------------

### High number of ID
ctcrw_lf_le_hID<- SDE$new(formulas = formulas,data = data_lf_le,type = "CTCRW",
                             response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                             other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
ctcrw_lf_le_hID$update_lambda(new_lambda)

#fit
ctcrw_lf_le_hID$fit(method="BFGS")

### Low number of ID
ctcrw_lf_le_lID<- SDE$new(formulas = formulas,data = data_lf_le[data_lf_le$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
ctcrw_lf_le_lID$update_lambda(new_lambda)

#fit
ctcrw_lf_le_lID$fit(method="BFGS")



## High frequency high error -------------------

### High number of ID
ctcrw_hf_he_hID<- SDE$new(formulas = formulas,data = data_hf_he,type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
ctcrw_hf_he_hID$update_lambda(new_lambda)

#fit
ctcrw_hf_he_hID$fit(method="BFGS")

### Low number of ID
ctcrw_hf_he_lID<- SDE$new(formulas = formulas,data = data_hf_le[data_hf_he$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
new_lambda=c(exp(4),exp(4))
ctcrw_hf_he_lID$update_lambda(new_lambda)

#fit
ctcrw_hf_he_lID$fit(method="BFGS")


## High frequency low error -------------------

### High number of ID
ctcrw_hf_le_hID<- SDE$new(formulas = formulas,data = data_hf_le,type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
ctcrw_hf_le_hID$update_lambda(new_lambda)

#fit
ctcrw_hf_le_hID$fit(method="BFGS")

### Low number of ID
ctcrw_hf_le_lID<- SDE$new(formulas = formulas,data = data_hf_le[data_hf_le$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize re variances
new_lambda=c(exp(4),exp(4))
ctcrw_hf_le_lID$update_lambda(new_lambda)

#fit
ctcrw_hf_le_lID$fit(method="BFGS")




# Estimate from simulated data with CRCVM  --------------------------------------------


## Low frequency high error ----------------
formulas <- list(mu1=~1,mu2=~1,tau=~s(ID,bs="re"),
                 nu=~s(ID,bs="re"),omega=~te(theta,Ishore,k=SP_DF,bs="cs"))

### High number of ID

crcvm_lf_he_hID<- SDE$new(formulas = formulas,data = data_lf_he,type = "RACVM_SSM",
                response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
new_lambda=c(1/SIGMA_TAU_0^2,1/SIGMA_NU_0^2,lambda_splines)
crcvm_lf_he_hID$update_lambda(new_lambda)
new_map=list("log_lambda"=factor(c(1,2,NA,NA)))
crcvm_lf_he_hID$update_map(new_map)

#fit
crcvm_lf_he_hID$fit(method="BFGS")

### Low number of ID

crcvm_lf_he_lID<- SDE$new(formulas = formulas,data = data_lf_he[data_lf_he$ID %in% 1:N_ID_LOW,],type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_lf_he_lID$update_lambda(new_lambda)
crcvm_lf_he_lID$update_map(new_map)

#fit
crcvm_lf_he_lID$fit(method="BFGS")


## Low frequency low error -----------

### High number of ID

crcvm_lf_le_hID<- SDE$new(formulas = formulas,data = data_lf_le,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_lf_le_hID$update_lambda(new_lambda)
crcvm_lf_le_hID$update_map(new_map)

#fit
crcvm_lf_le_hID$fit(method="BFGS")

### Low number of ID

crcvm_lf_le_lID<- SDE$new(formulas = formulas,data = data_lf_le[data_lf_le$ID %in% 1:N_ID_LOW,],type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_lf_le_lID$update_lambda(new_lambda)
crcvm_lf_le_lID$update_map(new_map)

#fit
crcvm_lf_le_lID$fit(method="BFGS")

## High frequency high error --------

### High number of ID

crcvm_hf_he_hID<- SDE$new(formulas = formulas,data = data_hf_he,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_hf_he_hID$update_lambda(new_lambda)
crcvm_hf_he_hID$update_map(new_map)

#fit
crcvm_hf_he_hID$fit(method="BFGS")

### Low number of ID

crcvm_hf_he_lID<- SDE$new(formulas = formulas,data = data_hf_he[data_hf_he$ID %in% 1:N_ID_LOW,],type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_hf_he_lID$update_lambda(new_lambda)
crcvm_hf_he_lID$update_map(new_map)

#fit
crcvm_hf_he_lID$fit(method="BFGS")

## High frequency low error -------


crcvm_hf_le_hID<- SDE$new(formulas = formulas,data = data_hf_le,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_hf_le_hID$update_lambda(new_lambda)
crcvm_hf_le_hID$update_map(new_map)

#fit
crcvm_hf_le_hID$fit(method="BFGS")

### Low number of ID

crcvm_hf_le_lID<- SDE$new(formulas = formulas,data = data_hf_le[data_hf_le$ID %in% 1:N_ID_LOW,],type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_hf_le_lID$update_lambda(new_lambda)
crcvm_hf_le_lID$update_map(new_map)

#fit
crcvm_hf_le_lID$fit(method="BFGS")


# Write estimated parameters in csv files -----------------------------

write_estimates_csv=function(model,name) {
  
  #coeffs
  coeffs=rbind(model$coeff_re(),model$coeff_fe()) #re and fe coeffs
  
  coeff_names=rownames(coeffs)
  coeff_values=as.numeric(coeffs)
  
  #standard deviation of random effects
  sdev=model$sdev()
  sdev_names=rownames(sdev)
  sdev_values=as.numeric(sdev)
  
  #measurement error
  log_sigma_obs_value=as.list(model$tmb_rep(),what="Est")$log_sigma_obs
  
  #create dataframe
  coeffs_df=data.frame("coeff_name"=factor(c(coeff_names,sdev_names,"log_sigma_obs")),
                                       "estimate"=c(coeff_values,sdev_values,log_sigma_obs_value))
  
  write.csv(coeffs_df,name, row.names=FALSE)
  
}
## CTCRW -----------

hyper_params_file_name=sapply(strsplit(hyperparams_file,"\\."),'[',1)

### Low frequency high error ---

#### High number of ID
write_estimates_csv(ctcrw_lf_he_hID,paste("result_",domain_name,"ctcrw_lf_he_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

#### Low number of ID
write_estimates_csv(ctcrw_lf_he_lID,paste("result_",domain_name,"ctcrw_lf_he_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

### Low frequency low error ---

#### High number of ID
write_estimates_csv(ctcrw_lf_le_hID,paste("result_",domain_name,"ctcrw_lf_le_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))
#### Low number of ID
write_estimates_csv(ctcrw_lf_le_lID,paste("result_",domain_name,"ctcrw_lf_le_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))


### High frequency high error ----

#### High number of ID
write_estimates_csv(ctcrw_hf_he_hID,paste("result_",domain_name,"ctcrw_hf_he_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

#### Low number of ID
write_estimates_csv(ctcrw_hf_he_lID,paste("result_",domain_name,"ctcrw_hf_he_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

### High frequency low error ----

#### High number of ID
write_estimates_csv(ctcrw_hf_le_hID,paste("result_",domain_name,"ctcrw_hf_le_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))
#### Low number of ID
write_estimates_csv(ctcrw_hf_le_lID,paste("result_",domain_name,"ctcrw_hf_le_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))


# CRCVM


### Low frequency high error ---

#### High number of ID
write_estimates_csv(crcvm_lf_he_hID,paste("result_",domain_name,"crcvm_lf_he_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

#### Low number of ID
write_estimates_csv(crcvm_lf_he_lID,paste("result_",domain_name,"crcvm_lf_he_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

### Low frequency low error ---

#### High number of ID
write_estimates_csv(crcvm_lf_le_hID,paste("result_",domain_name,"crcvm_lf_le_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))
#### Low number of ID
write_estimates_csv(crcvm_lf_le_lID,paste("result_",domain_name,"crcvm_lf_le_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))


### High frequency high error ----

#### High number of ID
write_estimates_csv(crcvm_hf_he_hID,paste("result_",domain_name,"crcvm_hf_he_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

#### Low number of ID
write_estimates_csv(crcvm_hf_he_lID,paste("result_",domain_name,"crcvm_hf_he_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))

### High frequency low error ----

#### High number of ID
write_estimates_csv(crcvm_hf_le_hID,paste("result_",domain_name,"crcvm_hf_le_hID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))
#### Low number of ID
write_estimates_csv(crcvm_hf_le_lID,paste("result_",domain_name,"crcvm_hf_le_lID_","seed",seed,
                                          hyper_params_file_name,".csv",sep=""))



# Save spline surface estimates ---------------

write_surface=function(model,name) {
  
  
  plots=model$get_all_plots(link=list("Ishore"=(\(x) 1/x)))
  
  surface=plots$fe_omega_theta_Ishore
  
  Dshore_values = as.numeric(surface$x$data[[1]]$x)
  theta_values = as.numeric(surface$x$data[[1]]$y)
  omega_values = as.numeric(surface$x$data[[1]]$z)
  
  data_surface=expand.grid(Dshore = Dshore_values, theta = theta_values)
  data_surface=as.vector(omega_values)
  
  write.csv(data_surface,name, row.names=FALSE)
  
}
## Low frequency high error


write_surface(crcvm_lf_he_hID,"surface_crcvm_lf_he_hID.csv")
write_surface(crcvm_lf_he_lID,"surface_crcvm_lf_he_lID.csv")


# Low frequency low error

write_surface(crcvm_lf_le_hID,"surface_crcvm_lf_le_hID.csv")
write_surface(crcvm_lf_le_lID,"surface_crcvm_lf_le_lID.csv")

## High frequency high error
write_surface(crcvm_hf_he_hID,"surface_crcvm_hf_he_hID.csv")
write_surface(crcvm_hf_he_lID,"surface_crcvm_hf_he_lID.csv")

## High frequency low error
write_surface(crcvm_hf_le_hID,"surface_crcvm_hf_le_hID.csv")
write_surface(crcvm_hf_le_lID,"surface_crcvm_hf_le_lID.csv")
