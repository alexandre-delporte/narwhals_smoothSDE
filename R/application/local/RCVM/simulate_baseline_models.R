# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-30
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/application/CRCVM/simulate_baseline_models.R
#
# Script Description:
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console


library(sf)
library(ggplot2)
library(smoothSDE)

# FJORDS POLYGONS


# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
greenland_data_path <- file.path(par_dir,"Data", "Greenland")  

#get the coastline geometry from the geojson file
land<-st_read(file.path(greenland_data_path,"land_utm.shp"))



#RECTANGULAR DOMAIN

#rectangular domain 
large_rectangle_coords <- matrix(c(0, 0, 0, 10000, 10000, 10000, 10000, 0, 0, 0), ncol = 2, byrow = TRUE)
small_rectangle_coords <- matrix(c(1000,1000, 1000, 9000, 9000, 9000, 9000, 1000, 1000, 1000), ncol = 2, byrow = TRUE)

# Create an sf object representing the rectangles
large_rectangle <- st_sfc(st_polygon(list(large_rectangle_coords)))
large_rectangle<-st_sf(geometry = large_rectangle)

small_rectangle <- st_sfc(st_polygon(list(small_rectangle_coords)))
small_rectangle<-st_sf(geometry = small_rectangle)

#change its crs
rect<- st_sym_difference(large_rectangle, small_rectangle)




# INITIAL POSITIONS

n_samples=6

#generate random initial points in the rectangle
z0_rect=matrix(rep(NA,n_samples*2),ncol=2)
colnames(z0_rect)=c("x1","x2")

for (i in 1:n_samples) {
  #choose location uniformly in the domain
  x=c(runif(1,min=1.5,max=8.5),runif(1,min=1.5,max=8.5))
  z0_rect[i,]=x
}


#generate random initial points
N_ID=6

z0=matrix(rep(NA,N_ID*2),ncol=2)
colnames(z0)=c("x1","x2")

i=1
while (i<=N_ID) {
  #choose location uniformly in the map
  x=c(runif(1,min=440,max=515),runif(1,min=7780,max=7880))
  p=nearest_shore_point(st_point(x),land)
  Dshore=(x[1]-p[1])^2+(x[2]-p[2])^2
  #keep it as initial position if it is at least 50 metres away from the shore
  if (Dshore>2) {
    z0[i,]=x
    i=i+1
  }
}



#FUNCTIONS TO COMPUTE COVARIATES ALONG THE WAY

fDshore=function(z,v,p) {
  
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #distance to shore
  Dshore=sqrt(normal[1]^2+normal[2]^2)
  
  return (Dshore)
}

fExpShore=function(z,v,p) {
  
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #distance to shore
  Dshore=sqrt(normal[1]^2+normal[2]^2)
  if (Dshore>3){
    return (0)
  }
  else if (Dshore<0.05) {
    return (20)
  }
  else {
    return (1/Dshore)
  }
}

fangle=function(z,v,p) {
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #return angle
  return (signed_angle(normal,v))
}


atw=list("ExpShore"=fExpShore,"AngleNormal"=fangle)


# TIME STEPS
Tmax=10
delta=1/60/6
times=seq(0,Tmax,by=delta)
n=length(times)
data_reg=data.frame("AngleNormal"=rep(0,n*length(unique(dataBE2$ID))),
                    "ExpShore"=rep(1,n*length(unique(dataBE2$ID))),times=seq(0,Tmax,by=delta),ID=rep(unique(dataBE2$ID),each=n))




# SIMULATE WITH CONSTRAINTS
sample=baseline3$simulate(z0=z0,data=data_reg,atw=atw,land=land,omega_times = 1)


#plots

p=ggplot()+geom_sf(data=land$geometry,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample,aes(x,y,col=ID),size=0.6)+theme_minimal()

p_rect=ggplot()+geom_sf(data=rect,fill="grey")+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample_rect,aes(x,y,col=ID),size=0.6)+theme_minimal()


saveWidget(ggplotly(p1), file ="baseline3/sample_baseline3.html")
saveWidget(ggplotly(p1_rect), file ="baseline3/sample_baseline3_rect.html")
