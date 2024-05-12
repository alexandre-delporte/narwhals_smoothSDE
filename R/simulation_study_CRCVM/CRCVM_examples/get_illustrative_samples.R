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
set.seed(2)

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

#for spline estimation
library(mgcv)

par_dir=dirname(getwd())
source(file.path(par_dir,"CVM_functions.R"))  #get functions to simulate trajectories

# Fjords domain --------------

# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
greenland_data_path <- file.path(par_dir,"Data", "Greenland")

#get the land and coastline geometries from the geojson file
land_border<-st_read(file.path(greenland_data_path,"land_utm.shp"))


# Initial velocity and location in fjords -------------

v0=c(0,0)

#generate random initial points
x0=matrix(rep(NA,2),ncol=2)
colnames(x0)=c("x1","x2")

i=1
while (i<=1) {
  #choose location uniformly in the map
  x=c(runif(1,min=430,max=500),runif(1,min=7760,max=7900))
  p=nearest_shore_point(st_point(x),land_border)
  Dshore=(x[1]-p[1])^2+(x[2]-p[2])^2
  #keep it as initial position if it is at least 50 metres away from the shore
  if (Dshore>2) {
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
rect_border<- st_sym_difference(large_rectangle, small_rectangle)


# Circular domain around the initial position -------

# Create circle geometries
large_circle <- st_buffer(st_sfc(st_point(x0)), dist = 10)
small_circle <- st_buffer(st_sfc(st_point(x0)), dist = 8)

# Convert circles to polygons
large_circle_polygon <- st_cast(large_circle, "POLYGON")
small_circle_polygon <- st_cast(small_circle, "POLYGON")

# Create sf objects
large_circle_sf <- st_sf(geometry = large_circle_polygon)
small_circle_sf <- st_sf(geometry = small_circle_polygon)

# Calculate the symmetric difference to get the border
circ_border <- st_sym_difference(large_circle_sf, small_circle_sf)


#  Parametric smooth parameters of the CRCVM -----------

#constant tau
ftau_constant=function(cov_data,tau=1) {
  return (tau)
}

ftau_bump=function(cov_data,epsilon=pi/8,alpha=pi/8,tau0=0.5,tau1=4) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  tau=(tau0+(tau1-tau0)/2*(tanh((theta+(pi/2+epsilon))/alpha)-tanh((theta+(pi/2-epsilon))/alpha))+
         (tau1-tau0)/2*(tanh((theta+(-pi/2+epsilon))/alpha)-tanh((theta+(-pi/2-epsilon))/alpha)))
  
  return (tau)
}

ftau_bump_modified=function(cov_data,epsilon=pi/5,alpha=pi/8,D0=3,kappa=1,tau0=0.5,tau1=3) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  coeff=exp(-kappa*(Dshore/D0)^2)
  tau=(tau0+(tau1-tau0)/2*(tanh((theta+(pi/2+epsilon))/alpha)-tanh((theta+(pi/2-epsilon))/alpha))*coeff+
         (tau1-tau0)/2*(tanh((theta+(-pi/2+epsilon))/alpha)-tanh((theta+(-pi/2-epsilon))/alpha))*coeff)
  
  return (tau)
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

fomega_fast=function(cov_data,D0=0.05,omega0=90*pi,lambda=2,kappa=0.2) {
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
fnu_constant=function(cov_data,nu=4) {
  return (nu)
}

# Splines smooth for omega -------------------- 
n <- 100000
SP_DF=3

#sample theta and DistanceShore points in the domain
theta <- runif(n,-pi,pi);
DistanceShore <- runif(n,0.1,8);
samples=data.frame(theta=theta,DistanceShore=DistanceShore)

#define grid of values
theta_v <- seq(-pi,pi,length=30)
Dshore_v<- seq(0.1,8,length=30)
pr <- data.frame(theta=rep(theta_v,30),DistanceShore=rep(Dshore_v,rep(30,30)))


#points on the surface perturbed by gaussian noise
f <- fomega(samples)
y <- f+0.1*rnorm(n)


#fit with bivariate splines te
m1 <- gam(y~te(theta,DistanceShore,k=SP_DF),knots=list(theta=seq(-pi,pi,len=SP_DF),
                                                       DistanceShore=seq(0.1,8,len=SP_DF)))

ExpShore=1/DistanceShore

#fit with bivariate splines te
m2 <- gam(y~te(theta,ExpShore,k=SP_DF),knots=list(theta=seq(-pi,pi,len=SP_DF),
                                                  ExpShore=seq(1/8,1/0.1,len=SP_DF)))


fomega_splines=function(cov_data) {
  
  if (is.null(cov_data$ExpShore)){
    omega=predict(m1,newdata=cov_data)
  }
  else if (is.null(cov_data$DistanceShore)){
    omega=predict(m2,newdata=cov_data)
  }
  
  return(as.numeric(omega))
}

# Plot the different smooth functions ---------------------------------------

Dshore=seq(from=0.05,to=3,length.out=30)
theta=seq(from=-pi,to=pi,length.out=30)
grid<- as.data.frame(expand.grid(Dshore,theta))
colnames(grid)=c("DistanceShore","theta")
tau_bump_modified=matrix(ftau_bump_modified(grid),30,30)
omega=matrix(fomega(grid),30,30)
omega_fast=matrix(fomega_fast(grid),30,30)
omega_splines=matrix(fomega_splines(grid),30,30)
tau_bump=ftau_bump(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(0,100)))
omega_close=fomega(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(0.1,100)))
omega_far=fomega(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(1,100)))
omega_splines_close=fomega_splines(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(0.1,100)))
omega_splines_far=fomega_splines(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(1,100)))

#perspective plots
# Open a PNG graphics device
png(file.path("smooth_omega1.png"))
persp(Dshore,theta,omega,axes=T,xlab="Distance to shore",ylab="Theta",zlab="Omega")
dev.off()

png(file.path("smooth_omega2.png"))
persp(Dshore,theta,omega,theta=270,axes=T,xlab="Distance to shore",ylab="Theta",zlab="Omega")
dev.off()

png(file.path("smooth_splines_omega1.png"))
persp(Dshore,theta,omega_splines,axes=T,xlab="Distance to shore",ylab="Theta",zlab="Omega")
dev.off()

png(file.path("smooth_splines_omega2.png"))
persp(Dshore,theta,omega_splines,theta=270,axes=T,xlab="Distance to shore",ylab="Theta",zlab="Omega")
dev.off()


plot_tau_bump=ggplot()+geom_line(aes(x=seq(from=-pi,to=pi,length.out=100),y=tau_bump))+ylab("Tau")+xlab("Theta")


plot_tau_bump_modified <- plot_ly(type = "surface",
                         contours = list(z = list(show = TRUE, start = min(tau_bump_modified), end = max(tau_bump_modified), size = 5,color="black")),
                         x = ~theta,y = ~Dshore,z=~tau_bump_modified,
                         colors=colorRamp(c("blue", "lightblue", "chartreuse3", "yellow", "red")))%>%
  layout(title=paste("Tau"),
         scene=list(xaxis = list(title = "theta",showgrid = F),
                    yaxis = list(title = "Distance to shore",showgrid = F),zaxis = list(title = "tau"))) 

plot_omega <- plot_ly(type = "surface",
                           contours = list(z = list(show = TRUE, start = min(omega), end = max(omega), size = 5,color="black")),
                           x = ~theta,y = ~Dshore,z=~omega,
                           colors=colorRamp(c("blue", "lightblue", "chartreuse3", "yellow", "red")))%>%
  layout(title=paste("Omega"),
         scene=list(xaxis = list(title = "Theta",showgrid = F),
                    yaxis = list(title = "Distance to shore",showgrid = F),zaxis = list(title = "omega"))) 

plot_omega_close=ggplot()+geom_line(aes(x=seq(from=-pi,to=pi,length.out=100),y=omega_close))+ylab("Omega")+xlab("Theta")

plot_omega_far=ggplot()+geom_line(aes(x=seq(from=-pi,to=pi,length.out=100),y=omega_far))+ylab("Omega")+xlab("Theta")


plot_omega_splines_close=ggplot()+geom_line(aes(x=seq(from=-pi,to=pi,length.out=100),y=omega_splines_close))+ylab("Omega")+xlab("Theta")

plot_omega_splines_far=ggplot()+geom_line(aes(x=seq(from=-pi,to=pi,length.out=100),y=omega_splines_far))+ylab("Omega")+xlab("Theta")

saveWidget(plot_omega, file ="smooth_omega_standard.html")
saveWidget(plot_tau_bump_modified, file ="smooth_tau_bump_modified.html")
ggsave("smooth_tau_bump.png",plot=plot_tau_bump,width=10,height=5)
ggsave("smooth_omega_close.png",plot=plot_omega_close,width=10,height=5)
ggsave("smooth_omega_far.png",plot=plot_omega_far,width=10,height=5)
ggsave("smooth_omega_splines_close.png",plot=plot_omega_splines_close,width=10,height=5)
ggsave("smooth_omega_splines_far.png",plot=plot_omega_splines_far,width=10,height=5)

# Simulate contrained CRCVM in fjords -------------------

# standard CRCVM in fjords
res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                             times=seq(0,5*24,by=1/60),land=land_border,verbose=FALSE)

data_sim_fjords_standard=res$sim
data_sim_fjords_standard$ID=factor(rep(1,length(data_sim_fjords_standard$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_fjords_standard=cbind(data_sim_fjords_standard,data_shore)

# spline standard CRCVM in fjords
res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,5*24,by=1/60),land=land_border,verbose=FALSE)

data_sim_fjords_splines=res$sim
data_sim_fjords_splines$ID=factor(rep(1,length(data_sim_fjords_splines$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_fjords_splines=cbind(data_sim_fjords_splines,data_shore)


#CRCVM persistent along the shore in fjords
res=sim_theta_CRCVM(ftau=ftau_bump,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                         times=seq(0,5*24,by=1/60),land=land_border,verbose=FALSE)
data_sim_fjords_pers=res$sim
data_sim_fjords_pers$ID=factor(rep(1,length(data_sim_fjords_pers$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_fjords_pers=cbind(data_sim_fjords_pers,data_shore)

#subsamble to 5 minutes
data_sim_fjords_standard=data_sim_fjords_standard[seq(1,length(data_sim_fjords_standard$time),by=5),]
data_sim_fjords_splines=data_sim_fjords_splines[seq(1,length(data_sim_fjords_splines$time),by=5),]
data_sim_fjords_pers=data_sim_fjords_pers[seq(1,length(data_sim_fjords_pers$time),by=5),]

# Simulate constrained CRCVM in rectangle ----------------

# standard CRCVM in rectangle
res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,5*24,by=1/60),land=rect_border,verbose=FALSE)

data_sim_rect_standard=res$sim
data_sim_rect_standard$ID=factor(rep(1,length(data_sim_rect_standard$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_rect_standard=cbind(data_sim_rect_standard,data_shore)


# standard spline CRCVM in rectangle
res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,5*24,by=1/60),land=rect_border,verbose=FALSE)

data_sim_rect_splines=res$sim
data_sim_rect_splines$ID=factor(rep(1,length(data_sim_rect_splines$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_rect_splines=cbind(data_sim_rect_splines,data_shore)

# #CRCVM persistent along the boundary in rectangle
res=sim_theta_CRCVM(ftau=ftau_bump,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,5*24,by=1/60),land=rect_border,verbose=FALSE)

data_sim_rect_pers=res$sim
data_sim_rect_pers$ID=factor(rep(1,length(data_sim_rect_pers$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_rect_pers=cbind(data_sim_rect_pers,data_shore)

data_sim_rect_standard=data_sim_rect_standard[seq(1,length(data_sim_rect_standard$time),by=5),]
data_sim_rect_splines=data_sim_rect_splines[seq(1,length(data_sim_rect_splines$time),by=5),]
data_sim_rect_pers=data_sim_rect_pers[seq(1,length(data_sim_rect_pers$time),by=5),]

# Simulate constrained CRCVM in circle ----------------

res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,5*24,by=1/60),land=circ_border,verbose=FALSE)

data_sim_circ_standard=res$sim
data_sim_circ_standard$ID=factor(rep(1,length(data_sim_circ_standard$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_circ_standard=cbind(data_sim_circ_standard,data_shore)

res=sim_theta_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,5*24,by=1/60),land=circ_border,verbose=FALSE)

data_sim_circ_splines=res$sim
data_sim_circ_splines$ID=factor(rep(1,length(data_sim_circ_splines$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_circ_splines=cbind(data_sim_circ_splines,data_shore)

res=sim_theta_CRCVM(ftau=ftau_bump,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,5*24,by=1/60),land=circ_border,verbose=FALSE)

data_sim_circ_pers=res$sim
data_sim_circ_pers$ID=factor(rep(1,length(data_sim_circ_pers$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_circ_pers=cbind(data_sim_circ_pers,data_shore)


data_sim_circ_standard=data_sim_circ_standard[seq(1,length(data_sim_circ_standard$time),by=5),]
data_sim_circ_splines=data_sim_circ_splines[seq(1,length(data_sim_circ_splines$time),by=5),]
data_sim_circ_pers=data_sim_circ_pers[seq(1,length(data_sim_circ_pers$time),by=5),]

# Simulate standard RCVM ---------------------

data_sim=sim_RACVM(mu1=0,mu2=0,beta=1,sigma=2*4/sqrt(pi*1),omega=0.5,v0=v0,x0=x0[1,],times=seq(0,5*24,by=5/60),log_sigma_obs=NULL,verbose=FALSE,keep_true_pos=FALSE) 
data_sim$ID=factor(rep(1,length(data_sim$y1)))



#Zoom on Fjords for the plot --------------------

x.min= min(data_sim_fjords_pers$y1)-1
x.max = max(data_sim_fjords_pers$y1)+1
y.min = min(data_sim_fjords_pers$y2)-1
y.max = max(data_sim_fjords_pers$y2)+1

bbox <- st_bbox(c(xmin=x.min, xmax = x.max, ymin =y.min, ymax = y.max), 
                crs = st_crs(land_border))

bbox_polygon <- st_as_sfc(st_bbox(bbox)
                          , crs = st_crs(land_border))

cropped_land_pers<- land_border%>% st_intersection(bbox_polygon)

x.min= min(data_sim_fjords_standard$y1)-1
x.max = max(data_sim_fjords_standard$y1)+1
y.min = min(data_sim_fjords_standard$y2)-1
y.max = max(data_sim_fjords_standard$y2)+1

bbox <- st_bbox(c(xmin=x.min, xmax = x.max, ymin =y.min, ymax = y.max), 
                crs = st_crs(land_border))

bbox_polygon <- st_as_sfc(st_bbox(bbox)
                          , crs = st_crs(land_border))

cropped_land_standard<- land_border%>% st_intersection(bbox_polygon)

x.min= min(data_sim_fjords_splines$y1)-1
x.max = max(data_sim_fjords_splines$y1)+1
y.min = min(data_sim_fjords_splines$y2)-1
y.max = max(data_sim_fjords_splines$y2)+1

bbox <- st_bbox(c(xmin=x.min, xmax = x.max, ymin =y.min, ymax = y.max), 
                crs = st_crs(land_border))

bbox_polygon <- st_as_sfc(st_bbox(bbox)
                          , crs = st_crs(land_border))

cropped_land_splines<- land_border%>% st_intersection(bbox_polygon)



# Create the plots ----------------------
plot_illust_fjords_standard=ggplot()+geom_sf(data=cropped_land_standard$geometry,fill="grey")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_path(data=data_sim_fjords_standard,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_point(data=data_sim_fjords_standard,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_fjords_standard%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_fjords_splines=ggplot()+geom_sf(data=cropped_land_splines$geometry,fill="grey")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_path(data=data_sim_fjords_splines,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_point(data=data_sim_fjords_splines,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_fjords_splines%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_fjords_pers=ggplot()+geom_sf(data=cropped_land_pers$geometry,fill="grey")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_path(data=data_sim_fjords_pers,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_point(data=data_sim_fjords_pers,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_fjords_pers%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")


plot_illust_standard=ggplot()+
  geom_path(data=data_sim,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_point(data=data_sim,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")


plot_illust_rect_standard=ggplot()+geom_sf(data=rect_border$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=data_sim_rect_standard,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_path(data=data_sim_rect_standard,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_rect_standard%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_rect_splines=ggplot()+geom_sf(data=rect_border$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=data_sim_rect_splines,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_path(data=data_sim_rect_splines,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_rect_splines%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_rect_pers=ggplot()+geom_sf(data=rect_border$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=data_sim_rect_pers,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_path(data=data_sim_rect_pers,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_rect_pers%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_circ_standard=ggplot()+geom_sf(data=circ_border$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=data_sim_circ_standard,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_path(data=data_sim_circ_standard,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_circ_standard%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_circ_splines=ggplot()+geom_sf(data=circ_border$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=data_sim_circ_splines,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_path(data=data_sim_circ_splines,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_circ_splines%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

plot_illust_circ_pers=ggplot()+geom_sf(data=circ_border$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=data_sim_circ_pers,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_path(data=data_sim_circ_pers,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim_circ_pers%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")


#SAVE

ggsave(filename="illustrative_sample_fjords_standard.png",plot=plot_illust_fjords_standard,width=10,height=5)
ggsave(filename="illustrative_sample_fjords_splines.png",plot=plot_illust_fjords_splines,width=10,height=5)
ggsave(filename="illustrative_sample_fjords_pers.png",plot=plot_illust_fjords_pers,width=10,height=5)
ggsave(filename="illustrative_sample_standard.png",plot=plot_illust_standard,width=10,height=5)
ggsave(filename="illustrative_sample_rect_standard.png",plot=plot_illust_rect_standard,width=10,height=5)
ggsave(filename="illustrative_sample_rect_splines.png",plot=plot_illust_rect_splines,width=10,height=5)
ggsave(filename="illustrative_sample_rect_pers.png",plot=plot_illust_rect_pers,width=10,height=5)
ggsave(filename="illustrative_sample_circ_standard.png",plot=plot_illust_circ_standard,width=10,height=5)
ggsave(filename="illustrative_sample_circ_splines.png",plot=plot_illust_circ_splines,width=10,height=5)
ggsave(filename="illustrative_sample_circ_pers.png",plot=plot_illust_circ_pers,width=10,height=5)

