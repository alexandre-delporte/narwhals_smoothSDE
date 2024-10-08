#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

domain_name="rect"
par_dir=dirname(getwd())
set_up_file=paste("set_up_",domain_name,".R",sep="")
source(file.path(par_dir,set_up_file))     #get set up for simulation study
source(file.path(dirname(par_dir),"CVM_functions.R"))  #get functions to simulate trajectories
library(smoothSDE)          #sde models
library(foreach)            #foreach loop
library(doParallel)         #parallel computing

seed= 14
set.seed(seed)



# Generate samples ---------------

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


data=foreach (i=1:N_ID,.combine='rbind',.packages=c("progress","MASS","sf","mgcv")) %dopar% {
  
  set.seed((seed-1)*N_ID+i)
  
  #constant nu
  fnu_constant=function(cov_data) {
    return (exp(true_log_nu[i]))
  }
  
  ftau_constant=function(cov_data) {
    return (exp(true_log_tau[i]))
  }
  res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,
                      log_sigma_obs=log(SIGMA_OBS),v0=v0,x0=x0[i,],times=times,land=border,verbose=FALSE)
  
  data_sim=res$sim
  data_sim$ID=factor(rep(i,length(data_sim$y1)))
  data_shore=res$shore[,c("p1","p2")]
  cbind(data_sim,data_shore)
}

#stop cluster
stopCluster(cl)





# Points that reached land ---------------
count=0
for (id in unique(data$ID)) {
  sub_data=data[data$ID==id,]
  if (nrow(sub_data) < n_obs) {
    count=count+1
    cat("ID",id,"reached land","\n",sep=" ")
  }
}

cat(count/N_ID*100,"percent of the samples reached land")

#subsample
data=data[seq(1,n_obs*N_ID,by=BY),]


# Save plot of the trajectories ------------

plot=ggplot()+geom_sf(data=border$geometry,fill="grey")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_point(data=data,mapping=aes(y1,y2,col=ID),size=0.1)+
  geom_path(data=data,mapping=aes(y1,y2,color=ID),size=0.1)+
  geom_point(data = data%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

ggsave(filename=paste("plot_",domain_name,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km",seed,".png",sep=""),plot=plot,width=10,height=5)

if (count>0) {
  stop("Stop : at least one trajectory reached the shore.")
}

# Add covariates to simulation -----------

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
  #data has at least columns "time","y1" "y2", "p1","p2"
  #this function new data frame with three more columns "theta", "DistanceShore" and "ExpShore"
  
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
    
    #matrix of normal vectors
    normal=as.matrix(sub_data[2:n_sub,c("y1","y2")]-sub_data[2:n_sub,c("p1","p2")])
    
    #angle between velocity and normal vector
    theta_coast=signed_angle(normal,vexp_df)
    
    #adjust lengths 
    theta_coast=c(theta_coast,1) 
    
    new_data[sub_ind,"theta"]=theta_coast
  }
  
  #distance to shore
  Dshore=sqrt((new_data$y1-new_data$p1)^2+(new_data$y2-new_data$p2)^2)
  new_data$DistanceShore=Dshore
  new_data$ExpShore=1/Dshore
  
  return(new_data)
}

data=add_covs(data)

#only estimate if observed Distances to shore include the interval [0.5,5]

if (min(data$DistanceShore)>D_low | max(data$DistanceShore)<D_up) {
  stop("Observed distance to shore are irrelevant with the fixed knots")
}

# Estimate from simulated data ------------------

formulas <- list(mu1=~1,mu2=~1,tau=~s(ID,bs="re"),
                 nu=~s(ID,bs="re"),omega=~te(theta,DistanceShore,k=SP_DF,bs="cs"))

knots=list("omega"=knots_DistanceShore)

crcvm<- SDE$new(formulas = formulas,data = data,type = "RACVM",
                    response = c("y1","y2"),par0 = PAR0,fixpar=c("mu1","mu2"),
                    other_data=list("H"=H),knots=knots)

#initialize smoothing penalties and re variances
new_lambda=c(1,1,1,1)
crcvm$update_lambda(new_lambda)
#fit
crcvm$fit(method="BFGS")

# Get estimated coefficients --------------


coeffs=rbind(crcvm$coeff_re(),crcvm$coeff_fe()) #re and fe coeffs

coeff_names=rownames(coeffs)
coeff_values=as.numeric(coeffs)

#standard deviation of random effects
sdev=crcvm$sdev()
sdev_names=rownames(sdev)
sdev_values=as.numeric(sdev)

#bind dataframe
coeffs_df=data.frame("coeff_name"=factor(c(coeff_names,sdev_names)),"estimate"=c(coeff_values,sdev_values))

#true values of the coeffs
true_df=data.frame("true"=c(tau_re,nu_re,sp_coeff_Dshore[2:SP_DF[1]*SP_DF[2]],0,0,
                            log(TAU_0),log(NU_0),sp_coeff_Dshore[1],SIGMA_TAU,SIGMA_NU,1/sqrt(m1$sp)))

coeffs_df=cbind(coeffs_df,true_df)

write.csv(coeffs_df,paste("result_",domain_name,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_DistanceShore",seed,".csv",sep=""), row.names=FALSE)

