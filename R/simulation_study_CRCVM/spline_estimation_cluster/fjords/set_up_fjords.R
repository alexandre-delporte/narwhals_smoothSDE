# HEADER ----------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-03
#
# Script Name:    set_up_fjords.R
#
# Script Description: Set up of the simulation study for the CRCVM within 
# Scoresby Sound fjords system
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
hyperparams_file <- "hyperparams_set1.txt"
path=here("R","simulation_study_CRCVM","spline_estimation_cluster","fjords")

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
  if (!(is_in_land(st_point(x),border))) {
    p=nearest_boundary_point(st_point(x),border)
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

fomega_cubic=function(cov_data,a,D0,D1,sigma_theta,sigma_D,b){
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  if (is.null(Dshore)){
    Dshore=1/cov_data$Ishore
  }
  
  a*theta*(theta-pi/2)*(theta+pi/2)*exp(-Dshore/D0)/Dshore+
    b*(exp(-1/2*(((theta+pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2))-
         exp(-1/2*(((theta-pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2)))
}


fomega=function(cov_data) {fomega_cubic(cov_data,A,D0,D1,SIGMA_THETA,SIGMA_D,B)} 


# Approximation of smooth omega with tensor splines -----------------------

n <- 100000                           #number of observations

theta <- runif(n,-pi,pi)            #sample theta
DistanceShore <- runif(n,D_LOW,D_UP)     #sample DistanceShore

samples=data.frame(theta=theta,DistanceShore=DistanceShore)  

#define grid of values
theta_v <- seq(-pi,pi,length=30)
Dshore_v<- seq(0.1,7,length=30)
pr <- data.frame(theta=rep(theta_v,30),DistanceShore=rep(Dshore_v,rep(30,30)))

# true values of the function over this grid
truth <- matrix(fomega(pr),30,30)

#points on the surface perturbed by gaussian noise
f <- fomega(samples)
y <- f+0.1*rnorm(n)

#plot true function
persp(theta_v,Dshore_v,truth);
title("truth")


Ishore=ifelse(DistanceShore<D_LOW,1/D_LOW,ifelse(DistanceShore>D_UP,0,1/DistanceShore))

#fit with bivariate splines te of Ishore
m <- gam(y~te(theta,Ishore,k=SP_DF,bs="cs"))


#visualize
vis.gam(m);title("tensor product")


#get splines coefficients=
sp_coeff=m$coefficients


#get smooothing penalties coefficients
lambda_splines=m$sp


#define spline smooth parameter omega
fomega_splines=function(cov_data) {
  
    omega=predict(m,newdata=cov_data)
  
  return(as.numeric(omega))
}

#k-th basis spline function for omega
fomega_spline_basis=function(cov_data,k=5) {
  
  Xp=predict(m2, newdata = cov_data, type = "lpmatrix")
  coeffs=rep(0,ncol(Xp))
  coeffs[k]=1
  return(as.numeric(Xp%*%coeffs))
}

