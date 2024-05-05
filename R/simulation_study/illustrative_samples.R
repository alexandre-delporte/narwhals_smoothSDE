# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-04
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/simulation_study/illustrative_samples.R
#
# Script Description: illustrative samples of the SDE models for the paper.
# 
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

#seed for reproducibility 
set.seed(42)

#parallel computing
library(foreach)
library(doParallel)

#for pipeline operator
library(dplyr)

#plots
library(ggplot2)
library(plotly)
library(htmlwidgets)

#data managing
library(tidyr)

#progress bar
library(progress)

#for land polygons
library(sf)

par_dir=dirname(getwd())
source(file.path(par_dir,"CVM_functions.R"))  #get functions to simulate trajectories

# Fjords domain --------------

# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
greenland_data_path <- file.path(par_dir,"Data", "Greenland")

#get the land and coastline geometries from the geojson file
land<-st_read(file.path(greenland_data_path,"land_utm.shp"))


# Initial velocity and location in fjords -------------

v0=c(0,0)

#generate random initial points
x0=matrix(rep(NA,2),ncol=2)
colnames(x0)=c("x1","x2")

i=1
while (i<=1) {
  #choose location uniformly in the map
  x=c(runif(1,min=430,max=500),runif(1,min=7760,max=7900))
  p=nearest_shore_point(st_point(x),land)
  Dshore=(x[1]-p[1])^2+(x[2]-p[2])^2
  #keep it as initial position if it is at least 50 metres away from the shore
  if (Dshore>0.5) {
    x0[i,]=x
    i=i+1
  }
}


# Rectangular domain around the initial location ---------------------
large_rectangle_coords <- matrix(c(x0[1,1]-10, x0[1,2]-10, x0[1,1]-10, x0[1,2]+10, x0[1,1]+10, x0[1,2]+10, x0[1,1]+10, x0[1,2]-10, x0[1,1]-10,  x0[1,2]-10), ncol = 2, byrow = TRUE)
small_rectangle_coords <-  matrix(c(x0[1,1]-9, x0[1,2]-9, x0[1,1]-9, x0[1,2]+9, x0[1,1]+9, x0[1,2]+9, x0[1,1]+9, x0[1,2]-9, x0[1,1]-9,x0[1,2]-9), ncol = 2, byrow = TRUE)

# Create an sf object representing the rectangles
large_rectangle <- st_sfc(st_polygon(list(large_rectangle_coords)))
large_rectangle<-st_sf(geometry = large_rectangle)

small_rectangle <- st_sfc(st_polygon(list(small_rectangle_coords)))
small_rectangle<-st_sf(geometry = small_rectangle)

#change its crs
domain<- st_sym_difference(large_rectangle, small_rectangle)



#  Smooth parameters of the RACVM and CRCVM -----------

#constant tau
ftau_constant=function(cov_data,tau=1,sigma_re=0.1) {
  b=rnorm(1,mean=0,sd=sigma_re)
  return (tau+b)
}

#constant omega
fomega_constant=function(cov_data,omega=0) {
  return (omega)
}

fomega=function(cov_data,D0=0.3,omega0=60*pi/2,lambda=2,kappa=0.2) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  if (is.null(Dshore)){
    Dshore=1/cov_data$ExpShore
  }
  coeff=exp(-kappa*(Dshore/D0)^2)
  omega=omega0/2*(tanh(lambda*(theta+pi/2-atanh(-0.9)/lambda))+tanh(lambda*(theta-(pi/2-atanh(-0.9)/lambda))))*coeff
  return(omega)
}

#constant nu
fnu_constant=function(cov_data,nu=4,sigma_re=0.3) {
  b=rnorm(1,mean=0,sd=sigma_re)
  return (nu+b)
}


# Simulate contrained CRCVM in fjords -------------------

res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],times=seq(0,5*24,by=1/60),land=land,verbose=FALSE)

data_sim_fjords=res$sim
data_sim_fjords$ID=factor(rep(1,length(data_sim_fjords$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_fjords=cbind(data_sim_fjords,data_shore)

#subsamble to 5 minutes
data_sim_fjords=data_sim_fjords[seq(1,length(data_sim_fjords$time),by=5),]

# Simulate constrained CRCVM in rectangle ----------------

res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],times=seq(0,5*24,by=1/60),land=domain,verbose=FALSE)

data_sim_rect=res$sim
data_sim_rect$ID=factor(rep(1,length(data_sim_rect$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_rect=cbind(data_sim_rect,data_shore)

data_sim_rect=data_sim_rect[seq(1,length(data_sim_rect$time),by=5),]

# Simulate standard RCVM ---------------------

data_sim=sim_RACVM(mu1=0,mu2=0,beta=1,sigma=2*4/sqrt(pi*1),omega=0.5,v0=v0,x0=x0[1,],times=seq(0,5*24,by=1/60),log_sigma_obs=NULL,verbose=FALSE,keep_true_pos=FALSE) 
data_sim$ID=factor(rep(1,length(data_sim$y1)))

data_sim=data_sim[seq(1,length(data_sim$time),by=5),]


#Zoom on Fjords for the plot --------------------

x.min= min(data_sim_fjords$y1)-1
x.max = max(data_sim_fjords$y1)+1
y.min = min(data_sim_fjords$y2)-1
y.max = max(data_sim_fjords$y2)+1

bbox <- st_bbox(c(xmin=x.min, xmax = x.max, ymin =y.min, ymax = y.max), 
                crs = st_crs(land))

bbox_polygon <- st_as_sfc(st_bbox(bbox)
                          , crs = st_crs(land))

cropped_land<- land%>% st_intersection(bbox_polygon)



# Create the plots ----------------------
plot_illust_fjords=ggplot()+geom_sf(data=cropped_land$geometry,fill="grey")+coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_path(data=data_sim_fjords,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_point(data=data_sim_fjords,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_fjords%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_standard=ggplot()+
  geom_path(data=data_sim,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_point(data=data_sim,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")


plot_illust_rect=ggplot()+geom_sf(data=domain$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=data_sim_rect,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_path(data=data_sim_rect,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_rect%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

#SAVE

ggsave(filename="illustrative_sample_fjords.png",plot=plot_illust_fjords,width=10,height=5)
ggsave(filename="illustrative_sample_standard.png",plot=plot_illust_standard,width=10,height=5)
ggsave(filename="illustrative_sample_rect.png",plot=plot_illust_rect,width=10,height=5)

