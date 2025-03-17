# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2025 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2025-01-16
#
# Script Name:    
#
# Script Description:
#
#
# SETUP ------------------------------------


set.seed(42)                #seed for reproductibility


library(ggplot2)            #plots
library(plotly)             #html plots
library(sf)                 #geometry for constraints
library(mgcv)               #tensor splines
library(smoothSDE)          #to compute nearest shore point


# READ HYPERPARAMETERS FILE

# Read hyperparameters from a user specified file
hyperparams_file="hyperparams_set2_24.txt"
path=here("R","simulation_study","cluster","CTCRW","fjords")

# Read the file and evaluate each line
lines <- readLines(file.path(path,hyperparams_file))
for (line in lines) {
  eval(parse(text = line))
}

# Land polygons --------------------


# Set the path to the directory containing the data
greenland_data_path <- here("Data","preprocessed_data","greenland")

#get the land and coastline geometries from the geojson file
border<-st_read(file.path(greenland_data_path,"updated_scoresby_sound_utm.shp"))
border<- st_transform(border,crs="+init=EPSG:32626 +units=km")

v0=c(0,0)

#generate random initial points
x0=matrix(rep(NA,N_ID_HIGH*2),ncol=2)
colnames(x0)=c("x1","x2")


i=1
while (i<=N_ID_HIGH) {
  #choose location uniformly in the map
  x=c(runif(1,min=430,max=500),runif(1,min=7760,max=7900))
  if (!(is_in_border(st_point(x),border))) {
    p=nearest_boundary_points(matrix(x,ncol=2),border)
    Dshore=(x[1]-p[1])^2+(x[2]-p[2])^2
    #keep it as initial position if it is at least 50 metres away from the shore
    if (Dshore>DMIN && Dshore<DMAX) {
      x0[i,]=x
      i=i+1
    }
  }
}

# Time steps ----------------


#high frequency time points
times=seq(0,TMAX,by=DELTA)

n_obs=length(times)-1

# Definition of parameters tau and nu ----------------
tau_re=rnorm(N_ID_HIGH,mean=0,sd=SIGMA_TAU)
nu_re=rnorm(N_ID_HIGH,mean=0,sd=SIGMA_NU)
true_log_tau=tau_re+log(TAU_0)
true_log_nu=nu_re+log(NU_0)


# Defintion of smooth parameter omega ----------


fomega_cubic=function(cov_data,boundary_distance_name="BoundaryDistance",
                      boundary_angle_name="BoundaryAngle",
                      a,b,D0,D1,sigma_D,sigma_theta){
  Dshore=cov_data[[boundary_distance_name]]
  theta=cov_data[[boundary_angle_name]]
  
  a*theta*(theta-pi/2)*(theta+pi/2)*exp(-Dshore/D0)/Dshore+
    b*(exp(-1/2*(((theta+pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2))-
         exp(-1/2*(((theta-pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2)))
}


fomega=function(cov_data,boundary_distance_name="BoundaryDistance",
                boundary_angle_name="BoundaryAngle") {
  fomega_cubic(cov_data,boundary_distance_name,
               boundary_angle_name,
               A,B,D0,D1,SIGMA_D,SIGMA_THETA)
}

