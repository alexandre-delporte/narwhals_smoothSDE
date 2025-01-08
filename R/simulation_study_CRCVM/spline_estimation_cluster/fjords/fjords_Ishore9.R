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


seed= 9
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
  res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega_cubic,fnu=fnu_constant,
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

cat(count/N_ID*100,"percent of the samples reached land")

if (count>0) {
  stop("At least one trajectory reached land.")
}


# add noise
low_noise=rmvn(n_obs,rep(0,2),diag(rep(SIGMA_OBS_LOW^2,2)))
high_noise=rmvn(n_obs,rep(0,2),diag(rep(SIGMA_OBS_HIGH^2,2)))
observed_data_le=data_sim
observed_data_le[,c("y1","y2")]=data_sim[,c("y1","y2")]+low_noise
observed_data_he=data_sim
observed_data_he[,c("y1","y2")]=data_sim[,c("y1","y2")]+high_noise

# subsample
data_lf_he=observed_data_he[seq(1,length(data_sim$time),by=BY_LF),]
data_lf_le=observed_data_le[seq(1,length(data_sim$time),by=BY_LF),]
data_hf_le=observed_data_le[seq(1,length(data_sim$time),by=BY_HF),]
data_hf_he=observed_data_he[seq(1,length(data_sim$time),by=BY_HF),]


# Save plot of the trajectories

plot=ggplot()+geom_sf(data=border$geometry,fill="lightgrey",border="black")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_point(data=data,mapping=aes(y1,y2,col=ID),size=0.1)+
  geom_path(data=data,mapping=aes(y1,y2,color=ID),size=0.1)+
  geom_point(data = data%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

ggsave(filename=paste("simulated_trajectories_","seed",seed,sapply(strsplit(hyperparams_file,"\\.")),".png",sep=""),plot=plot,width=10,height=5)

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
    new_data[sub_ind,"Ishore"]=ifelse(Dshore>D_LOW,1/Dshore,1/D_LOW)
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
                other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

#initialize re variances
new_lambda=c(1/SIGMA_TAU_0^2,1/SIGMA_NU_0^2)
ctcrw_lf_he_hID$update_lambda(new_lambda)

#fit
ctcrw_lf_he_hID$fit(method="BFGS")

### Low number of ID
ctcrw_lf_he_lID<- SDE$new(formulas = formulas,data = data_lf_he[data_lf_he$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                             response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                             other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

#initialize re variances
ctcrw_lf_he_lID$update_lambda(new_lambda)

#fit
ctcrw_lf_he_lID$fit(method="BFGS")

## Low frequency low error ------------------

### High number of ID
ctcrw_lf_le_hID<- SDE$new(formulas = formulas,data = data_lf_le,type = "CTCRW",
                             response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                             other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

#initialize re variances
ctcrw_lf_le_hID$update_lambda(new_lambda)

#fit
ctcrw_lf_le_hID$fit(method="BFGS")

### Low number of ID
ctcrw_lf_le_lID<- SDE$new(formulas = formulas,data = data_lf_le[data_lf_le$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

#initialize re variances
ctcrw_lf_le_lID$update_lambda(new_lambda)

#fit
ctcrw_lf_le_lID$fit(method="BFGS")



## High frequency high error -------------------

### High number of ID
ctcrw_hf_he_hID<- SDE$new(formulas = formulas,data = data_hf_he,type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

#initialize re variances
ctcrw_hf_he_hID$update_lambda(new_lambda)

#fit
ctcrw_hf_he_hID$fit(method="BFGS")

### Low number of ID
ctcrw_hf_he_lID<- SDE$new(formulas = formulas,data = data_hf_le[data_hf_he$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

#initialize re variances
new_lambda=c(exp(4),exp(4))
ctcrw_hf_he_lID$update_lambda(new_lambda)

#fit
ctcrw_hf_he_lID$fit(method="BFGS")


## High frequency low error -------------------

### High number of ID
ctcrw_hf_le_hID<- SDE$new(formulas = formulas,data = data_hf_le,type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

#initialize re variances
ctcrw_hf_le_hID$update_lambda(new_lambda)

#fit
ctcrw_hf_le_hID$fit(method="BFGS")

### Low number of ID
ctcrw_hf_le_lID<- SDE$new(formulas = formulas,data = data_hf_le[data_hf_le$ID %in% 1:N_ID_LOW,],type = "CTCRW",
                          response = c("y1","y2"),par0 = PAR0[1:4],fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS_0)))

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

crcvm_lf_he_hID<- SDE$new(formulas = formulas,data = data_lf_he_hID,type = "RACVM_SSM",
                response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
new_lambda=c(1/SIGMA_TAU_0^2,1/SIGMA_NU_0^2,lambda_splines)
crcvm_lf_he_hID$update_lambda(new_lambda)

#fit
crcvm_lf_he_hID$fit(method="BFGS")

### Low number of ID

crcvm_lf_he_lID<- SDE$new(formulas = formulas,data = data_lf_he_lID,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_lf_he_lID$update_lambda(new_lambda)

#fit
crcvm_lf_he_lID$fit(method="BFGS")


## Low frequency low error -----------

### High number of ID

crcvm_lf_le_hID<- SDE$new(formulas = formulas,data = data_lf_le_hID,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_lf_le_hID$update_lambda(new_lambda)

#fit
crcvm_lf_le_hID$fit(method="BFGS")

### Low number of ID

crcvm_lf_le_lID<- SDE$new(formulas = formulas,data = data_lf_le_lID,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_lf_le_lID$update_lambda(new_lambda)

#fit
crcvm_lf_le_lID$fit(method="BFGS")

## High frequency high error --------

### High number of ID

crcvm_hf_he_hID<- SDE$new(formulas = formulas,data = data_hf_he_hID,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_hf_he_hID$update_lambda(new_lambda)

#fit
crcvm_hf_he_hID$fit(method="BFGS")

### Low number of ID

crcvm_hf_he_lID<- SDE$new(formulas = formulas,data = data_hf_he_lID,type = "RACVM_SSM",
                          response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                          other_data=list("log_sigma_obs0"=log(SIGMA_OBS0)))

#initialize smoothing penalties and re variances
crcvm_hf_he_lID$update_lambda(new_lambda)

#fit
crcvm_hf_he_lID$fit(method="BFGS")




# Write estimated parameters in csv files -----------------------------

write_estimates_csv=function(model,name) {
  
  coeffs=rbind(model$coeff_re(),model$coeff_fe()) #re and fe coeffs
  
  coeff_names=rownames(coeffs)
  coeff_values=as.numeric(coeffs)
  
  #standard deviation of random effects
  sdev=model$sdev()
  sdev_names=rownames(sdev_ctcrw_lf_he_hID)
  sdev_values=as.numeric(sdev)
  
  #create dataframe
  coeffs_df=data.frame("coeff_name"=factor(c(coeff_names,sdev_names)),
                                       "estimate"=c(coeff_values,sdev_values))
  
  write.csv(coeffs_df,name, row.names=FALSE)
  
}
## CTCRW -----------

### Low frequency high error ---

#### High number of ID
write_estimates_csv(ctcrw_lf_he_hID,paste("result_",domain_name,"ctcrw_lf_he_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

#### Low number of ID
write_estimates_csv(ctcrw_lf_he_lID,paste("result_",domain_name,"ctcrw_lf_he_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

### Low frequency low error ---

#### High number of ID
write_estimates_csv(ctcrw_lf_le_hID,paste("result_",domain_name,"ctcrw_lf_le_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))
#### Low number of ID
write_estimates_csv(ctcrw_lf_le_lID,paste("result_",domain_name,"ctcrw_lf_le_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))


### High frequency high error ----

#### High number of ID
write_estimates_csv(ctcrw_hf_he_hID,paste("result_",domain_name,"ctcrw_hf_he_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

#### Low number of ID
write_estimates_csv(ctcrw_hf_he_lID,paste("result_",domain_name,"ctcrw_hf_he_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

### High frequency low error ----

#### High number of ID
write_estimates_csv(ctcrw_hf_le_hID,paste("result_",domain_name,"ctcrw_hf_le_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))
#### Low number of ID
write_estimates_csv(ctcrw_hf_le_lID,paste("result_",domain_name,"ctcrw_hf_le_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))


# CRCVM


### Low frequency high error ---

#### High number of ID
write_estimates_csv(crcvm_lf_he_hID,paste("result_",domain_name,"crcvm_lf_he_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

#### Low number of ID
write_estimates_csv(crcvm_lf_he_lID,paste("result_",domain_name,"crcvm_lf_he_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

### Low frequency low error ---

#### High number of ID
write_estimates_csv(crcvm_lf_le_hID,paste("result_",domain_name,"crcvm_lf_le_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))
#### Low number of ID
write_estimates_csv(crcvm_lf_le_lID,paste("result_",domain_name,"crcvm_lf_le_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))


### High frequency high error ----

#### High number of ID
write_estimates_csv(crcvm_hf_he_hID,paste("result_",domain_name,"crcvm_hf_he_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

#### Low number of ID
write_estimates_csv(crcvm_hf_he_lID,paste("result_",domain_name,"crcvm_hf_he_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))

### High frequency low error ----

#### High number of ID
write_estimates_csv(crcvm_hf_le_hID,paste("result_",domain_name,"crcvm_hf_le_hID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))
#### Low number of ID
write_estimates_csv(crcvm_hf_le_lID,paste("result_",domain_name,"crcvm_hf_le_lID_","seed",seed,
                                          sapply(strsplit(hyperparams_file,"\\.")),".csv",sep=""))





# Get estimated parameters for RACVM-----------------------------


## Low frequency high error ------------

### High number of ID
coeffs_crcvm_lf_he_hID=rbind(crcvm_lf_he_hID$coeff_re(),crcvm_lf_he_hID$coeff_fe()) #re and fe coeffs

coeff_names_crcvm_hID=rownames(coeffs_lf_he_hID)
coeff_values_crcvm_lf_he_hID=as.numeric(coeffs_crcvm_lf_he_hID)

#standard deviation of random effects
sdev_crcvm_lf_he_hID=sdev_crcvm_lf_he_hID$sdev()
sdev_names_crcvm=rownames(sdev_crcvm_lf_he_hID)
sdev_values_crcvm_lf_he_hID=as.numeric(sdev_crcvm_lf_he_hID)

#create dataframe
coeffs_df_crcvm_lf_he_hID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_hID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_lf_he_hID,sdev_values_crcvm_lf_he_hID))

### Low number of ID

coeffs_crcvm_lf_he_lID=rbind(crcvm_lf_he_lID$coeff_re(),crcvm_lf_he_lID$coeff_fe()) #re and fe coeffs

coeff_names_crcvm_lID=rownames(coeffs_lf_he_lID)
coeff_values_crcvm_lf_he_lID=as.numeric(coeffs_crcvm_lf_he_lID)

#standard deviation of random effects
sdev_crcvm_lf_he_lID=sdev_crcvm_lf_he_lID$sdev()
sdev_values_crcvm_lf_he_lID=as.numeric(sdev_crcvm_lf_he_lID)

#create dataframe
coeffs_df_crcvm_lf_he_lID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_lID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_lf_he_lID,sdev_values_crcvm_lf_he_lID))



## Low frequency low error ------------

### High number of ID
coeffs_crcvm_lf_le_hID=rbind(crcvm_lf_le_hID$coeff_re(),crcvm_lf_le_hID$coeff_fe()) #re and fe coeffs

coeff_values_crcvm_lf_le_hID=as.numeric(coeffs_crcvm_lf_le_hID)

#standard deviation of random effects
sdev_crcvm_lf_le_hID=sdev_crcvm_lf_le_hID$sdev()
sdev_values_crcvm_lf_le_hID=as.numeric(sdev_crcvm_lf_le_hID)

#create dataframe
coeffs_df_crcvm_lf_le_hID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_hID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_lf_le_hID,sdev_values_crcvm_lf_le_hID))

### Low number of ID

coeffs_crcvm_lf_le_lID=rbind(crcvm_lf_le_lID$coeff_re(),crcvm_lf_le_lID$coeff_fe()) #re and fe coeffs

coeff_values_crcvm_lf_le_lID=as.numeric(coeffs_crcvm_lf_le_lID)

#standard deviation of random effects
sdev_crcvm_lf_le_lID=sdev_crcvm_lf_le_lID$sdev()
sdev_values_crcvm_lf_le_lID=as.numeric(sdev_crcvm_lf_le_lID)

#create dataframe
coeffs_df_crcvm_lf_le_lID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_lID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_lf_le_lID,sdev_values_crcvm_lf_le_lID))




## High frequency high error ------------

### High number of ID
coeffs_crcvm_hf_he_hID=rbind(crcvm_hf_he_hID$coeff_re(),crcvm_hf_he_hID$coeff_fe()) #re and fe coeffs

coeff_values_crcvm_hf_he_hID=as.numeric(coeffs_crcvm_hf_he_hID)

#standard deviation of random effects
sdev_crcvm_hf_he_hID=sdev_crcvm_hf_he_hID$sdev()
sdev_values_crcvm_hf_he_hID=as.numeric(sdev_crcvm_hf_he_hID)

#create dataframe
coeffs_df_crcvm_hf_he_hID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_hID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_hf_he_hID,sdev_values_crcvm_hf_he_hID))

### Low number of ID

coeffs_crcvm_hf_he_lID=rbind(crcvm_hf_he_lID$coeff_re(),crcvm_hf_he_lID$coeff_fe()) #re and fe coeffs

coeff_values_crcvm_hf_he_lID=as.numeric(coeffs_crcvm_hf_he_lID)

#standard deviation of random effects
sdev_crcvm_hf_he_lID=sdev_crcvm_hf_he_lID$sdev()
sdev_values_crcvm_hf_he_lID=as.numeric(sdev_crcvm_hf_he_lID)

#create dataframe
coeffs_df_crcvm_hf_he_lID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_lID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_hf_he_lID,sdev_values_crcvm_hf_he_lID))



## High frequency low error ------------

### High number of ID
coeffs_crcvm_hf_le_hID=rbind(crcvm_hf_le_hID$coeff_re(),crcvm_hf_le_hID$coeff_fe()) #re and fe coeffs

coeff_values_crcvm_hf_le_hID=as.numeric(coeffs_crcvm_hf_le_hID)

#standard deviation of random effects
sdev_crcvm_hf_le_hID=sdev_crcvm_hf_le_hID$sdev()
sdev_values_crcvm_hf_le_hID=as.numeric(sdev_crcvm_hf_le_hID)

#create dataframe
coeffs_df_crcvm_hf_le_hID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_hID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_hf_le_hID,sdev_values_crcvm_hf_le_hID))

### Low number of ID

coeffs_crcvm_hf_le_lID=rbind(crcvm_hf_le_lID$coeff_re(),crcvm_hf_le_lID$coeff_fe()) #re and fe coeffs

coeff_values_crcvm_hf_le_lID=as.numeric(coeffs_crcvm_hf_le_lID)

#standard deviation of random effects
sdev_crcvm_hf_le_lID=sdev_crcvm_hf_le_lID$sdev()
sdev_values_crcvm_hf_le_lID=as.numeric(sdev_crcvm_hf_le_lID)

#create dataframe
coeffs_df_crcvm_hf_le_lID=data.frame("coeff_name"=factor(c(coeff_names_crcvm_lID,sdev_names_crcvm)),
                                     "estimate"=c(coeff_values_crcvm_hf_le_lID,sdev_values_crcvm_hf_le_lID))



# Plot spline estimates ---------------

## Low frequency high error

plots_crcvm_lf_he_lID=crcvm_lf_he_hID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))

plots_crcvm_lf_he_hID=crcvm_lf_he_lID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))


# Low frequency low error


plots_crcvm_lf_le_lID=crcvm_lf_le_hID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))

plots_crcvm_lf_le_hID=crcvm_lf_le_lID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))


## High frequency high error
plots_crcvm_hf_he_lID=crcvm_hf_he_lID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))

plots_crcvm_hf_he_hID=crcvm_hf_he_lID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))


## High frequency low error
plots_crcvm_hf_le_lID=crcvm_hf_le_lID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))

plots_crcvm_hf_le_hID=crcvm_hf_le_lID$plot_par(var="theta",par_name="omega",covs=data.frame("ID"="Asgeir","Ishore"=c(1/0.1,1/1,1/2.5)))

