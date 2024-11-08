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
# Script Description: illustrative samples of the CRCVM SDE models for the paper.
# 
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

#seed for reproducibility 
set.seed(99)

#parallel computing
library(foreach)
library(doParallel)

#for pipeline operator
library(dplyr)

#plots
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(ggpubr)

#data managing
library(tidyr)

#progress bar
library(progress)

#for land polygons
library(sf)

#for spline estimation
library(mgcv)

#to get root of git repo
library(here)

source(file.path(here("R","simulation_study_CRCVM","CVM_functions.R")))  #get functions to simulate trajectories

# Fjords domain --------------

# Set the path to the directory containing the data
greenland_data_path <- here("Data","preprocessed_data","greenland")

#get the land and coastline geometries from the geojson file
land_border<-st_read(file.path(greenland_data_path,"updated_scoresby_sound_utm.shp"))


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
small_circle <- st_buffer(st_sfc(st_point(x0)), dist = 9)

# Convert circles to polygons
large_circle_polygon <- st_cast(large_circle, "POLYGON")
small_circle_polygon <- st_cast(small_circle, "POLYGON")

# Create sf objects
large_circle_sf <- st_sf(geometry = large_circle_polygon)
small_circle_sf <- st_sf(geometry = small_circle_polygon)

# Calculate the symmetric difference to get the border
circ_border <- st_sym_difference(large_circle_sf, small_circle_sf)


# Create small circles borders within the big circle
circle1 <- st_buffer(st_sfc(st_point(x0-0.0005*x0)), dist = 2)

# Convert circles to polygons
circle_polygon1 <- st_cast(circle1, "POLYGON")


# Create sf objects
circle_sf1 <- st_sf(geometry = circle_polygon1)


# Create small circles borders within the big circle
circle2 <- st_buffer(st_sfc(st_point(x0+0.0005*x0)), dist = 2)

# Convert circles to polygons
circle_polygon2 <- st_cast(circle2, "POLYGON")


# Create sf objects
circle_sf2 <- st_sf(geometry = circle_polygon2)

combined_circ_border <- st_union(st_union(circ_border, circle_sf1), circle_sf2)




#  Parametric smooth parameters of the CRCVM -----------

#constant tau
ftau_constant=function(cov_data,tau=2) {
  return (tau)
}

ftau_bump=function(cov_data,epsilon=pi/12,alpha=pi/8,tau0=0.5,tau1=4) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  tau=(tau0+(tau1-tau0)/2*(tanh((theta+(pi/2+epsilon))/alpha)-tanh((theta+(pi/2-epsilon))/alpha))+
         (tau1-tau0)/2*(tanh((theta+(-pi/2+epsilon))/alpha)-tanh((theta+(-pi/2-epsilon))/alpha)))
  
  return (tau)
}

ftau_bump_rect=function(cov_data,epsilon=pi/10,alpha=pi/10,tau0=0.5,tau1=4) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  tau=(tau0+(tau1-tau0)/2*(tanh((theta+(pi/2+epsilon))/alpha)-tanh((theta+(pi/2-epsilon))/alpha))+
         (tau1-tau0)/2*(tanh((theta+(-pi/2+epsilon))/alpha)-tanh((theta+(-pi/2-epsilon))/alpha)))
  
  return (tau)
}

ftau_bump_modified=function(cov_data,epsilon=pi/8,alpha=pi/8,D0=3,kappa=1,tau0=0.5,tau1=4) {
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

fomega_fast=function(cov_data,D0=0.3,omega0=60*pi,lambda=2,kappa=0.2) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  if (is.null(Dshore)){
    Dshore=1/cov_data$ExpShore
  }
  coeff=exp(-kappa*(Dshore/D0)^2)
  omega=omega0/2*(tanh(lambda*(theta+pi/2-atanh(-0.9)/lambda))+tanh(lambda*(theta-(pi/2-atanh(-0.9)/lambda))))*coeff
  return(omega)
}

fomega_tortuous=function(cov_data,D0=0.3,omega0=60*pi/2,lambda=2,kappa=0.2) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  if (is.null(Dshore)){
    Dshore=1/cov_data$ExpShore
  }
  coeff=exp(-kappa*(Dshore/D0)^2)
  omega=omega0/2*(tanh(lambda*(theta+pi/2-atanh(-0.5)/lambda))+tanh(lambda*(theta-(pi/2-atanh(-0.5)/lambda))))*coeff
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
DistanceShore <- runif(n,0.2,8);
samples=data.frame(theta=theta,DistanceShore=DistanceShore)

#define grid of values
theta_v <- seq(-pi,pi,length=30)
Dshore_v<- seq(0.2,8,length=30)
pr <- data.frame(theta=rep(theta_v,30),DistanceShore=rep(Dshore_v,rep(30,30)))


#points on the surface perturbed by gaussian noise
f <- fomega(samples)
y <- f+0.1*rnorm(n)


#fit with bivariate splines te
m1 <- gam(y~te(theta,DistanceShore,k=SP_DF),knots=list(theta=seq(-pi,pi,len=SP_DF),
                                                       DistanceShore=seq(0.2,8,len=SP_DF)))

ExpShore=1/DistanceShore

#fit with bivariate splines te
m2 <- gam(y~te(theta,ExpShore,k=SP_DF),knots=list(theta=seq(-pi,pi,len=SP_DF),
                                                  ExpShore=seq(1/8,1/0.2,len=SP_DF)))


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

Dshore=seq(from=0.05,to=4,length.out=30)
theta=seq(from=-pi,to=pi,length.out=30)
grid<- as.data.frame(expand.grid(Dshore,theta))
colnames(grid)=c("DistanceShore","theta")
tau_bump_modified=matrix(ftau_bump_modified(grid),30,30)
omega=matrix(fomega(grid),30,30)
omega_fast=matrix(fomega_fast(grid),30,30)
omega_splines=matrix(fomega_splines(grid),30,30)
tau_bump=ftau_bump(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=500),"DistanceShore"=rep(0,500)))
omega_close=fomega(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(0.1,100)))
omega_far=fomega(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(1,100)))
omega_splines_close=fomega_splines(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(0.1,100)))
omega_splines_far=fomega_splines(cov_data=data.frame("theta"=seq(from=-pi,to=pi,length.out=100),"DistanceShore"=rep(1,100)))

#perspective plots
# Open a PNG graphics device
png(file.path("combined_persp_plots.png"), width = 800, height = 800)

# Set up the layout: 2 rows and 2 columns
par(mfrow = c(2, 2))

# First plot
persp(Dshore, theta, omega,phi=15, 
      axes = FALSE, # turn off default axes
      main = "Parametric function",
      cex.main = 2,ticktype="simple")  # Make the main title bigger
axis(1, labels = FALSE, tick = FALSE)
mtext("Distance to shore", side = 1, line = 1, cex = 1.5)
mtext(expression(omega), side = 2, line = 1, cex = 1.5, las = 2)  # Rotate y-axis label

# Second plot
persp(Dshore, theta, omega, phi=15,
      theta = -90, 
      axes = FALSE, 
      main = "Parametric function",
      cex.main = 2,ticktype="simple")  # Make the main title bigger
mtext(expression(Theta), side = 1, line = 1, cex = 1.5)
mtext(expression(omega), side = 2, line = 1, cex = 1.5, las = 2)  # Rotate y-axis label

# Third plot
persp(Dshore, theta, omega_splines, phi=15,
      axes = FALSE, 
      main = "Weighted splines",
      cex.main = 2,ticktype="simple")  # Make the main title bigger
mtext("Distance to shore", side = 1, line = 1, cex = 1.5)
mtext(expression(omega), side = 2, line = 1, cex = 1.5, las = 2)  # Rotate y-axis label

# Fourth plot
persp(Dshore, theta, omega_splines, phi=15,
      theta = -90, 
      axes = FALSE, 
      main = "Weighted splines",
      cex.main = 2,ticktype="simple")  # Make the main title bigger
mtext(expression(Theta), side = 1, line = 1, cex = 1.5)
mtext(expression(omega), side = 2, line = 1, cex = 1.5, las = 2)  # Rotate y-axis label

# Close the PNG device
dev.off()


plot_tau_bump <- ggplot() + 
  geom_line(aes(x = seq(from = -pi, to = pi, length.out = 500), y = tau_bump)) + 
  ylab(expression(tau)) + 
  xlab(expression(Theta)) + 
  geom_vline(xintercept = c(-pi/2, pi/2), linetype = "dashed", color = "red") + 
  annotate("text", x = -pi/2-0.3, y = 0.3, label = expression(-pi/2), color = "red", size = 6) + 
  annotate("text", x = pi/2+0.2, y = 0.3, label = expression(pi/2), color = "red", size = 6) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.title.y = element_text(angle = 0)  # Rotate the y-axis label and adjust horizontal justification
  )


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
ggsave("smooth_tau_bump.pdf",plot=plot_tau_bump,width=10,height=10)
ggsave("smooth_omega_close.png",plot=plot_omega_close,width=10,height=5)
ggsave("smooth_omega_far.png",plot=plot_omega_far,width=10,height=5)
ggsave("smooth_omega_splines_close.png",plot=plot_omega_splines_close,width=10,height=5)
ggsave("smooth_omega_splines_far.png",plot=plot_omega_splines_far,width=10,height=5)

# Simulate contrained CRCVM in fjords -------------------

TMAX=5*24
# standard CRCVM in fjords
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                             times=seq(0,TMAX,by=1/60),land=land_border,verbose=FALSE)

data_sim_fjords_standard=res$sim
data_sim_fjords_standard$ID=factor(rep(1,length(data_sim_fjords_standard$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_fjords_standard=cbind(data_sim_fjords_standard,data_shore)

# spline standard CRCVM in fjords
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=land_border,verbose=FALSE)

data_sim_fjords_splines=res$sim
data_sim_fjords_splines$ID=factor(rep(1,length(data_sim_fjords_splines$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_fjords_splines=cbind(data_sim_fjords_splines,data_shore)


#CRCVM persistent along the shore in fjords
set.seed(99)
res=sim_CRCVM(ftau=ftau_bump,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                         times=seq(0,TMAX,by=1/60),land=land_border,verbose=FALSE)
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
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega_fast,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=rect_border,verbose=FALSE)

data_sim_rect_standard=res$sim
data_sim_rect_standard$ID=factor(rep(1,length(data_sim_rect_standard$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_rect_standard=cbind(data_sim_rect_standard,data_shore)


# standard spline CRCVM in rectangle
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=rect_border,verbose=FALSE)

data_sim_rect_splines=res$sim
data_sim_rect_splines$ID=factor(rep(1,length(data_sim_rect_splines$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_rect_splines=cbind(data_sim_rect_splines,data_shore)

# #CRCVM persistent along the boundary in rectangle
set.seed(99)
res=sim_CRCVM(ftau=ftau_bump_rect,fomega=fomega_fast,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=rect_border,verbose=FALSE)

data_sim_rect_pers=res$sim
data_sim_rect_pers$ID=factor(rep(1,length(data_sim_rect_pers$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_rect_pers=cbind(data_sim_rect_pers,data_shore)

data_sim_rect_standard=data_sim_rect_standard[seq(1,length(data_sim_rect_standard$time),by=5),]
data_sim_rect_splines=data_sim_rect_splines[seq(1,length(data_sim_rect_splines$time),by=5),]
data_sim_rect_pers=data_sim_rect_pers[seq(1,length(data_sim_rect_pers$time),by=5),]

# Simulate constrained CRCVM in circle ----------------

#standard CRCVM in circle
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=circ_border,verbose=FALSE)

data_sim_circ_standard=res$sim
data_sim_circ_standard$ID=factor(rep(1,length(data_sim_circ_standard$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_circ_standard=cbind(data_sim_circ_standard,data_shore)

#spline CRCVM in circle
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=circ_border,verbose=FALSE)

data_sim_circ_splines=res$sim
data_sim_circ_splines$ID=factor(rep(1,length(data_sim_circ_splines$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_circ_splines=cbind(data_sim_circ_splines,data_shore)

#persistent CRCVM in circle
set.seed(99)
res=sim_CRCVM(ftau=ftau_bump,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=circ_border,verbose=FALSE)

data_sim_circ_pers=res$sim
data_sim_circ_pers$ID=factor(rep(1,length(data_sim_circ_pers$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_circ_pers=cbind(data_sim_circ_pers,data_shore)


data_sim_circ_standard=data_sim_circ_standard[seq(1,length(data_sim_circ_standard$time),by=5),]
data_sim_circ_splines=data_sim_circ_splines[seq(1,length(data_sim_circ_splines$time),by=5),]
data_sim_circ_pers=data_sim_circ_pers[seq(1,length(data_sim_circ_pers$time),by=5),]


# Simulate constrained CRCVM in combined circles ----------------

#standard CRCVM in combined circles
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=combined_circ_border,verbose=FALSE)

data_sim_combined_circ_standard=res$sim
data_sim_combined_circ_standard$ID=factor(rep(1,length(data_sim_combined_circ_standard$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_combined_circ_standard=cbind(data_sim_combined_circ_standard,data_shore)

#spline CRCVM in circle
set.seed(99)
res=sim_CRCVM(ftau=ftau_constant,fomega=fomega_splines,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=combined_circ_border,verbose=FALSE)

data_sim_combined_circ_splines=res$sim
data_sim_combined_circ_splines$ID=factor(rep(1,length(data_sim_combined_circ_splines$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_combined_circ_splines=cbind(data_sim_combined_circ_splines,data_shore)

#persistent CRCVM in circle
set.seed(99)
res=sim_CRCVM(ftau=ftau_bump,fomega=fomega,fnu=fnu_constant,v0=v0,x0=x0[1,],
                    times=seq(0,TMAX,by=1/60),land=combined_circ_border,verbose=FALSE)

data_sim_combined_circ_pers=res$sim
data_sim_combined_circ_pers$ID=factor(rep(1,length(data_sim_combined_circ_pers$y1)))
data_shore=res$shore[,c("p1","p2")]
data_sim_combined_circ_pers=cbind(data_sim_combined_circ_pers,data_shore)


data_sim_combined_circ_standard=data_sim_combined_circ_standard[seq(1,length(data_sim_combined_circ_standard$time),by=5),]
data_sim_combined_circ_splines=data_sim_combined_circ_splines[seq(1,length(data_sim_combined_circ_splines$time),by=5),]
data_sim_combined_circ_pers=data_sim_combined_circ_pers[seq(1,length(data_sim_combined_circ_pers$time),by=5),]


# Simulate standard RCVM ---------------------

data_sim=sim_RACVM(mu1=0,mu2=0,beta=1,sigma=2*4/sqrt(pi*1),omega=0.5,v0=v0,x0=x0[1,],times=seq(0,TMAX,by=5/60),log_sigma_obs=NULL,verbose=FALSE,keep_true_pos=FALSE) 
data_sim$ID=factor(rep(1,length(data_sim$y1)))



#Zoom on Fjords for the plot --------------------

x.min= min(c(data_sim_fjords_pers$y1,data_sim_fjords_standard$y1,data_sim_fjords_splines$y1))-1
x.max = max(c(data_sim_fjords_pers$y1,data_sim_fjords_standard$y1,data_sim_fjords_splines$y1))+1
y.min = min(c(data_sim_fjords_pers$y2,data_sim_fjords_standard$y2,data_sim_fjords_splines$y2))-1
y.max = max(c(data_sim_fjords_pers$y2,data_sim_fjords_standard$y2,data_sim_fjords_splines$y2))+1

bbox <- st_bbox(c(xmin=x.min, xmax = x.max, ymin =y.min, ymax = y.max), 
                crs = st_crs(land_border))

bbox_polygon <- st_as_sfc(st_bbox(bbox)
                          , crs = st_crs(land_border))

cropped_land<- land_border%>% st_intersection(bbox_polygon)



# Create the plots ----------------------

#standard unconstrained 
plot_illust_standard=ggplot()+
  geom_path(data=data_sim,mapping=aes(y1,y2,col=time),size=0.1)+
  geom_point(data=data_sim,mapping=aes(y1,y2,col=time),size=0.1)+
  scale_color_viridis_c(name="Time")+
  geom_point(data = data_sim%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")+theme_minimal()


# Ensure all geometries have the same CRS
st_crs(cropped_land$geometry) <- st_crs("+init=EPSG:32626 +units=km")
st_crs(rect_border$geometry) <- st_crs("+init=EPSG:32626 +units=km")
st_crs(circ_border$geometry) <- st_crs("+init=EPSG:32626 +units=km")

# Combine geometries into a single sf object
geometries <- bind_rows(
  st_as_sf(data.frame(geometry = cropped_land$geometry, domain = "Fjords")),
  st_as_sf(data.frame(geometry = rect_border$geometry, domain = "Rect")),
  st_as_sf(data.frame(geometry = circ_border$geometry, domain = "Circ"))
)

# Combine all data into one dataframe
data_combined <- bind_rows(
  data_sim_fjords_standard %>% mutate(domain = "Fjords", model = "Standard CRCVM"),
  data_sim_fjords_splines %>% mutate(domain = "Fjords", model = "Tortuous CRCVM"),
  data_sim_fjords_pers %>% mutate(domain = "Fjords", model = "Persistent CRCVM"),
  data_sim_rect_standard %>% mutate(domain = "Rect", model = "Standard CRCVM"),
  data_sim_rect_splines %>% mutate(domain = "Rect", model = "Tortuous CRCVM"),
  data_sim_rect_pers %>% mutate(domain = "Rect", model = "Persistent CRCVM"),
  data_sim_circ_standard %>% mutate(domain = "Circ", model = "Standard CRCVM"),
  data_sim_circ_splines %>% mutate(domain = "Circ", model = "Tortuous CRCVM"),
  data_sim_circ_pers %>% mutate(domain = "Circ", model = "Persistent CRCVM"),
  
)

# Function to create plots with individual scales
create_plot <- function(domain, model) {
  data_filtered <- data_combined %>% filter(domain == !!domain & model == !!model)
  
  ggplot() +
    geom_sf(data = geometries %>% filter(domain == !!domain), aes(geometry = geometry), fill = "grey", inherit.aes = FALSE) +
    coord_sf(datum = st_crs("+init=EPSG:32626 +units=km")) +
    geom_path(data = data_filtered, aes(x = y1, y = y2, col = time), size = 0.1) +
    geom_point(data = data_filtered, aes(x = y1, y = y2, col = time), size = 0.1) +
    geom_point(data = x0, aes(x = x1, y = x2), shape = 3, size = 2, col = "red") +
    scale_color_viridis_c(name = "Time") +
    theme_minimal() +
    xlab("x") + ylab("y") +
    theme(legend.position = "none")  # Remove legend from individual plots
}


# Get unique combinations of domain and model
unique_combinations <- expand.grid(model = unique(data_combined$model),domain = unique(data_combined$domain))

# Create a list of plots
plot_list <- lapply(seq_len(nrow(unique_combinations)), function(i) {
  domain <- unique_combinations$domain[i]
  model <- unique_combinations$model[i]
  create_plot(domain, model)
})

final_plot=ggarrange(plotlist=plot_list,
                     ncol=3, nrow=3, common.legend = TRUE, legend="bottom",
                     labels= c("Standard CRCVM", "Tortuous CRCVM", "Persistent CRCVM"))
# Print final combined plot
print(final_plot)

ggsave("final_combined_plot.png", plot = final_plot, width = 14, height = 14)
