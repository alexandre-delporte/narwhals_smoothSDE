# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2025 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2025-03-21
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/application/cluster/CRCVM/fit_baseline.R
#
# Script Description:
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space
.libPaths(c("~/R/library", .libPaths()))

library(smoothSDE)
library(ggplot2)
library(dplyr)
library(plotly)
library(mgcv)
library(here)
library(xtable)
library(sf)
library(doParallel)
library(foreach)
library(ggpubr)
library(fitdistrplus)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a hyperparameter file as an argument.")
}

hyperparams_file <- args[1]


# get name of hyperparameters file without extension
hyperparams_file_name=sapply(strsplit(hyperparams_file,"\\."),'[',1)

path=here("R","application","cluster","CRCVM")


# Read hyperparameters from a user specified file
lines <- readLines(file.path(path,hyperparams_file))
for (line in lines) {
  eval(parse(text = line))
}

#Auxiliary functions
source(file.path(path,"utility.R"))


# Set the path to the directory containing the data
par_dir=here()

narwhal_data_path <- file.path(par_dir,"Data","preprocessed_data","narwhals")  

#data before exposure with first 12h removed
dataBE12=read.csv(file.path(narwhal_data_path,"DataBE12.csv"), header = TRUE,dec = ".")

#data after exposure
dataAE=read.csv(file.path(narwhal_data_path,"DataAE.csv"), header = TRUE,dec = ".")



# Set the path to the directory containing the greenland data
greenland_data_path <- here("Data","preprocessed_data","greenland")


border200<-st_read(file.path(greenland_data_path,"shrunk200_scoresby_sound_utm.shp"))
border200 <- st_transform(border200, crs = "+init=EPSG:32626 +units=km")


#prepare data for interpolation
dataBE12_prep200=dataBE12[,c("ID","x","y","time","theta200","DistanceShore200")]

colnames(dataBE12_prep200)[colnames(dataBE12_prep200)=="DistanceShore200"]<-"BoundaryDistance"
colnames(dataBE12_prep200)[colnames(dataBE12_prep200)=="theta200"]<-"BoundaryAngle"

#remove points on land
dataBE12_prep200<-dataBE12_prep200[dataBE12_prep200$BoundaryDistance>0,]

n_step=N_STEP

#compute spline interpolation
interpolation_bas_200<-interpolate_BoundaryMetrics(dataBE12_prep200,response=c("x","y"),
                                                   border200,n_step=N_STEP,n_cores=6,df_ratio=DF_RATIO)

D0=quantile(dataBE12_prep200$BoundaryDistance,D0_QUANTILE)
D1=quantile(dataBE12_prep200$BoundaryDistance,D1_QUANTILE)
a=A*D0

baseline_crcvm200<- fit_baseline(dataBE12_prep200,interpolation_bas_200,
                                 b=B,a=a,sigma_D=SIGMA_D,D0=D0,D1=D1,sigma_theta=SIGMA_THETA,
                                 sigma_obs=SIGMA_OBS)




# Check baseline model

## Simulate trajectories

#initial narwhals positions (convert to matrix to avoid bugs)
z0_BE <- as.matrix(dataBE12_prep200[!duplicated(dataBE12_prep200$ID), c("x", "y")])

# Define functions to compute covariates along the way
fangle=function(z,v,p) {
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #return angle
  return (signed_angle(normal,v))
}


fDshore=function(z,v,p) {
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #distance to shore
  Dshore=sqrt(normal[1]^2+normal[2]^2)
  
  #fix a threshold to avoid error when we hit the boundary
  if (Dshore>0.005) {
    return(Dshore)
  }
  else {
    return(0.005)
  }
  
}

atw=list("BoundaryAngle"=fangle,"BoundaryDistance"=fDshore)


# Define time steps
delta=1/60/6
times_list=lapply(unique(dataBE12_prep200$ID),function(ID) {
  min=min(dataBE12[dataBE12_prep200$ID==ID,"time"])
  max=max(dataBE12[dataBE12_prep200$ID==ID,"time"])
  seq(min,max,by=delta)})
names(times_list)=unique(dataBE12_prep200$ID)
times_df <- data.frame(
  ID = rep(names(times_list), sapply(times_list, length)),
  time = unlist(times_list),row.names = NULL
)
n=length(times_df$time)
# Define covariates data (needed in the function simulate)
data_reg=data.frame("BoundaryDistance"=rep(1,n),
                    "BoundaryAngle"=rep(0,n))

data_reg=cbind(data_reg,times_df)


all_sim=list()

for (i in 1:N_SIM) {
  set.seed(i)
  baseline_crcvm200_trajectory=baseline_crcvm200$simulate(z0=z0_BE,
                                                          data=data_reg,atw=atw,land=border200,verbose=FALSE,n_cores=6)
  colnames(baseline_crcvm200_trajectory)[colnames(baseline_crcvm200_trajectory)=="z1"]<-"x"
  colnames(baseline_crcvm200_trajectory)[colnames(baseline_crcvm200_trajectory)=="z2"]<-"y"
  
  all_sim[[i]]<-baseline_crcvm200_trajectory
  
}



## Simulate time steps with log normal distribution
dtimes<-unlist(lapply(unique(dataBE12_prep200$ID),function(id) {
  sub_data<-dataBE12_prep200[dataBE12_prep200$ID==id,]
  dt=diff(sub_data$time)
  dt
}))

fit_lognorm <- fitdist(dtimes * 60, "lnorm")
fit_gamma <- fitdist(dtimes * 60, "gamma")
gofstat(list(fit_lognorm,fit_gamma))

# Extract estimated parameters
meanlog <- fit_lognorm$estimate["meanlog"]
sdlog <- fit_lognorm$estimate["sdlog"]

# Generate new time steps
dtimes_sim <- rlnorm(nrow(dataBE12_prep200)-6, meanlog, sdlog)

times_sim<-lapply(unique(dataBE12_prep200$ID),function(ID) {
  sub_time<- dataBE12_prep200[dataBE12_prep200$ID==ID,"time"]
  sub_n<-length(sub_time)
  min_time=min(sub_time)
  times<-rep(min_time,sub_n)+c(0,cumsum(dtimes_sim[1:(sub_n-1)]/60))
  return (times)})
names(times_sim)<-unique(dataBE12_prep200$ID)


all_sim_sub=list()

for (i in seq_along(all_sim)) {
  sim_data<-all_sim[[i]]
  interp_sim_data<-dataBE12_prep200[,c("time","ID","x","y")]
  interp_sim_data$time<-unlist(times_sim)
  
  for (id in unique(sim_data$ID)) {
    sub_sim_data<-sim_data[sim_data$ID==id,]
    sim_x_interp <- approx(sub_sim_data$time, sub_sim_data$x, xout = times_sim[[id]])$y
    sim_y_interp <- approx(sub_sim_data$time, sub_sim_data$y, xout = times_sim[[id]])$y
    sub_ind=interp_sim_data$ID==id
    interp_sim_data[sub_ind,"x"]<-sim_x_interp
    interp_sim_data[sub_ind,"y"]<-sim_y_interp
  }
  
  all_sim_sub[[i]]<-interp_sim_data
}

#subsample to 5 minutes and add measurement noise
all_sim_sub<-lapply(all_sim_sub,function(data) {
  
  data<-data[!(is.na(data$x)),]
  noise<-rmvn(nrow(data),rep(0,2),diag(rep(SIGMA_OBS^2,2)))
  
  data[,c("x","y")]=data[,c("x","y")]+noise
  return (data)
})


plot_list<-plot_trajectories(all_sim_sub,dataBE12_prep200)

lapply(plot_list,function(plot) {ggsave(plot = plot, paste0("sim",i,"_",hyperparams_file_name,".pdf"), 
                        width = 10, height = 6, units = "in", limitsize = FALSE)})rplus

all_boundary_data<- lapply(all_sim_sub,function(data) {
  get_BoundaryMetrics(data,
                      response=c("x","y"),
                      border=border200)
})

for (i in seq_along(all_sim_sub)) {
  
  all_sim_sub[[i]][,c("BoundaryDistance","BoundaryAngle")]<- all_boundary_data[[i]]
  all_sim_sub[[i]]$speed<- compute_speed(data=all_sim_sub[[i]])
}


dataBE12_prep200$speed<-compute_speed(dataBE12_prep200)

plots<-make_density_plots(all_sim_sub,dataBE12_prep200) 

final_densities_plots<-ggarrange(plotlist=list(plots$distance,plots$angle),
                                 ncol=2, nrow=1, common.legend = TRUE, legend="bottom", vjust=2)
ggsave(paste0("BoundaryDistance_density_",hyperparams_file_name,".pdf"),plot=plots$distance)
ggsave(paste0("BoundaryAngle_density_",hyperparams_file_name,".pdf"),plot=plots$angle)
ggsave(paste0("speed_density_",hyperparams_file_name,".pdf"),plot=plots$speed)




