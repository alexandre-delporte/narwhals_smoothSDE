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
cat("\014")                 
rm(list = ls())           
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

path=here("R","application","cluster","CRCVM","measurement_error_sensitivity")


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


make_final_baseline_tab(baseline_crcvm200,include_sigma_obs=FALSE,n_post=100,resp=TRUE)

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


all_sim_sub=list()

for (i in seq_along(all_sim)) {
  sim_data<-all_sim[[i]]
  interp_sim_data<-dataBE12_prep200[,c("time","ID","x","y")]
  
  for (id in unique(sim_data$ID)) {
    sub_sim_data<-sim_data[sim_data$ID==id,]
    sim_x_interp <- approx(sub_sim_data$time, sub_sim_data$x, 
                           xout = dataBE12_prep200[dataBE12_prep200$ID==id,"time"])$y
    sim_y_interp <- approx(sub_sim_data$time, sub_sim_data$y,
                           xout = dataBE12_prep200[dataBE12_prep200$ID==id,"time"])$y
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

percentage<-sum(unlist(lapply(all_sim_sub,function(data) {
  sum(tapply(data$BoundaryDistance == 0, data$ID, any))
})))/(N_SIM*6)
cat(percentage,"% of simulated trajectories reached boundary.","\n",sep=" ")

dataBE12_prep200$speed<-compute_speed(dataBE12_prep200)

plots<-make_density_plots(all_sim_sub,dataBE12_prep200) 

final_densities_plots<-ggarrange(plotlist=list(plots$distance,plots$angle),
                                 ncol=2, nrow=1, common.legend = TRUE, legend="bottom", vjust=2)
ggsave(paste0("BoundaryDistance_density_",hyperparams_file_name,".pdf"),plot=plots$distance)
ggsave(paste0("BoundaryAngle_density_",hyperparams_file_name,".pdf"),plot=plots$angle)
ggsave(paste0("speed_density_",hyperparams_file_name,".pdf"),plot=plots$speed)




results_200<-propagate_baseline_uncertainty(baseline_crcvm200,dataAE_prep200,
                                                   interpolation_resp_200,n_sim=100,
                                                   sigma_obs=SIGMA_OBS)
alpha_estimates<- results_200$estimates


index<- rowSums(alpha_estimates == 0) < 2
alpha_estimates_filtered<-alpha_estimates[index, ]
mean_alpha=list(alpha_tau=mean(alpha_estimates_filtered[,1]),
                alpha_nu=mean(alpha_estimates_filtered[,2]))
quant_alpha<- list(
  alpha_tau=quantile(alpha_estimates_filtered[,1],probs = c(0.025, 0.975)),
  alpha_nu=quantile(alpha_estimates_filtered[,2],probs = c(0.025, 0.975)))

parameter_names <- names(mean_alpha)
estimates <- mean_alpha
conf_intervals <- sapply(parameter_names, function(param) {
  sprintf("$[%.2f; %.2f]$", quant_alpha[[param]][1], 
          quant_alpha[[param]][2])
})


table_data <- data.frame(
  Parameter = 
    c("$\\alpha_{\\tau}$","$\\alpha_{\\nu}$"),
  Estimate = sprintf("$%.2f$", mean_alpha),
  CI = conf_intervals,
  stringsAsFactors = FALSE
)
xtab <- xtable(
  table_data,
  caption = "Response log linear estimations with propagated uncertainty",
  label = "table:response_log_linear_estimations_with_propagated_uncertainty"
)

print(
  xtab,type="latex",
  include.rownames = FALSE,
  sanitize.text.function = identity,  
  hline.after = c(-1, 0, nrow(table_data)),
  table.placement="H"
  
)



models_list_filtered<- results_200$models[index]
CI_plots=lapply(models_list_filtered,function(model) { 
  model$get_all_plots(baseline=NULL,xmin=xmin,xmax=xmax,
                      link=link,xlabel=xlabel,show_CI="pointwise")})


# Extract data for each estimated smooth parameter
tau_CI_data <- lapply(CI_plots, function(plot) {
  ggplot_build(plot$fe_tau_ExpShip)$data[[1]][,c("x","y")]
})
nu_CI_data <- lapply(CI_plots, function(plot) {
  ggplot_build(plot$fe_nu_ExpShip)$data[[1]][,c("x","y")]
})

# Get confidence intervals as quantiles of estimated parameters

# Combine the x columns from all data frames into a matrix
x_matrix_tau <- do.call(cbind, lapply(tau_CI_data, function(df) df$y))
x_matrix_nu <- do.call(cbind, lapply(nu_CI_data, function(df) df$y))
# Compute the 5% and 95% quantiles row-wise
quantiles_tau <- apply(x_matrix_tau, 1, 
                       function(row) quantile(row, probs = c(0.025, 0.975)))
quantiles_nu <- apply(x_matrix_nu, 1,
                      function(row) quantile(row, probs = c(0.025, 0.975)))

# Transpose and format the results into a data frame
quantiles_tau <- as.data.frame(t(quantiles_tau))
colnames(quantiles_tau) <- c("low", "up")
quantiles_nu <- as.data.frame(t(quantiles_nu))
colnames(quantiles_nu) <- c("low", "up")

# Add the "true" estimated smooth
plots=response_crcvm200$get_all_plots(
  baseline=baseline_crcvm200,xmin=xmin,xmax=xmax,
  link=link,xlabel=xlabel,show_CI="pointwise")
est_tau=ggplot_build(plots$fe_tau_ExpShip)$data[[1]][,c("x","y")]
est_nu=ggplot_build(plots$fe_nu_ExpShip)$data[[1]][,c("x","y")]
est_tau=cbind(est_tau,quantiles_tau)
est_nu=cbind(est_nu,quantiles_nu) 

est_tau$baseline=ggplot_build(plots$fe_tau_ExpShip)$data[[2]][,"y"]
est_nu$baseline=ggplot_build(plots$fe_nu_ExpShip)$data[[2]][,"y"]


plot_tau_response_crcvm200<- ggplot(data = est_tau) +
  geom_line(aes(x = x, y = y), color = "steelblue", size = 1) +
  geom_line(aes(x=x,y=baseline),color="steelblue",size=1,linetype="dashed")+
  geom_ribbon(aes(x = x, ymin = low, ymax = up), 
              fill = "lightblue", alpha = 0.5) +
  labs(
    x = "Distance to Ship (km)",
    y = expression(tau)) +
  theme_minimal() +theme(
    axis.title.x = element_text(size=22),
    axis.title.y = element_text(size=22, angle=0, vjust=0.5, margin = margin(r = 10)),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    legend.title = element_text(size=18),
    legend.text = element_text(size=18)
  )

# Enhanced plot for nu
plot_nu_response_crcvm200 <- ggplot(data = est_nu) +
  geom_line(aes(x = x, y = y), color = "darkorange", size = 1) +
  geom_line(aes(x,baseline),color="darkorange",size=1,linetype="dashed")+
  geom_ribbon(aes(x = x, ymin = low, ymax = up), 
              fill = "navajowhite", alpha = 0.5) +
  labs(
    x = "Distance to Ship (km)",
    y = expression(nu)) +
  theme_minimal() +theme(
    axis.title.x = element_text(size=22),
    axis.title.y = element_text(size=22, angle=0, vjust=0.5, margin = margin(r = 10)),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    legend.title = element_text(size=18),
    legend.text = element_text(size=18)
  )

ggsave(paste0("fe_tau_ExpShip_",hyperparams_file_name,".pdf"), plot=plot_tau_response_crcvm200,
       width=10, height=6, units="in")

ggsave(paste0("fe_tau_ExpShip_",hyperparams_file_name,".pdf"), plot=plot_nu_response_crcvm200, width=10, height=6, units="in")
plot_tau_response_crcvm200
plot_nu_response_crcvm200

make_recovery_distances_tab<-function(p) {
  
  parameter_names=names(mean_alpha)
  scaling=as.list(c(log(1-p),log(1+p)))
  names(scaling)=parameter_names
  estimates <- sapply(parameter_names,function(param) {
    mean_alpha[[param]]/scaling[[param]]})
  
  conf_intervals <- sapply(parameter_names, function(param) {
    sprintf("$[%.2f; %.2f]$", quant_alpha[[param]][1]/scaling[[param]], 
            quant_alpha[[param]][2]/scaling[[param]])
  })
  
  
  table_data <- data.frame(
    Parameter = 
      c("$D_{\\tau}^{ship}$","$D_{\\nu}^{ship}$"),
    Estimate = sprintf("$%.2f$", estimates),
    CI = conf_intervals,
    stringsAsFactors = FALSE
  )
  
  xtab <- xtable(
    table_data,
    caption = "Estimated recovery distances of the baseline values",
    label = paste0("table:recovery_distances_",p)
  )
  
  print(
    xtab,type="latex",
    include.rownames = FALSE,
    sanitize.text.function = identity,  
    hline.after = c(-1, 0, nrow(table_data)),
    table.placement="H"
    
  )
}

make_recovery_distances_tab(0.5)

make_recovery_distances_tab(0.1)

