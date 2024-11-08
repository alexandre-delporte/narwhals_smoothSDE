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
# Script Description: Set up of the simulation study for the CRCVM witin 
# a Scoresby Sound fjords system
#
#
# SETUP ------------------------------------

set.seed(42)                #seed for reproductibility


library(ggplot2)            #plots
library(plotly)             #html plots
library(sf)                 #geometry for constraints
library(mgcv)               #tensor splines




# READ HYPERPARAMETERS FILE

# Read hyperparameters from a user specified file
hyperparams_file <- "hyperparams_file.txt"  # Update this with a bash variable
path=here("R","simulation_study_CRCVM","spline_estimation")

# Read the file and evaluate each line
lines <- readLines(file.path(path,hyperparams_file))
for (line in lines) {
  eval(parse(text = line))
}

# Circular domain --------------------

# Create circle geometries
large_circle <- st_buffer(st_sfc(st_point(c(0,0))), dist = 10)
small_circle <- st_buffer(st_sfc(st_point(c(0,0))), dist = 8)

# Convert circles to polygons
large_circle_polygon <- st_cast(large_circle, "POLYGON")
small_circle_polygon <- st_cast(small_circle, "POLYGON")

# Create sf objects
large_circle_sf <- st_sf(geometry = large_circle_polygon)
small_circle_sf <- st_sf(geometry = small_circle_polygon)

# Calculate the symmetric difference to get the border
border <- st_sym_difference(large_circle_sf, small_circle_sf)


# Initial velocity and position -------------

v0=c(0,0)

#generate random initial points
x0=matrix(rep(NA,N_ID*2),ncol=2)
colnames(x0)=c("x1","x2")

for (i in 1:N_ID) {
  #choose location uniformly in the circle https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
  r=3*sqrt(runif(1))
  theta=runif(1)*2*pi
  x=c(r*cos(theta),r*sin(theta))
  x0[i,]=x
}


# Time steps ----------------


#high frequency time points
times=seq(0,TMAX,by=DELTA)

n_obs=length(times)-1

# Measurement error -------------

H=array(rep(SIGMA_OBS^2*diag(2),n_obs*N_ID/BY),dim=c(2,2,n_obs*N_ID/BY))


# Definition of parameters tau and nu ----------------
tau_re=rnorm(N_ID,mean=0,sd=SIGMA_TAU)
nu_re=rnorm(N_ID,mean=0,sd=SIGMA_NU)
true_log_tau=tau_re+log(TAU_0)
true_log_nu=nu_re+log(NU_0)



# Defintion of smooth parameter omega ----------

fomega=function(cov_data,D0=0.5,omega0=20*pi,lambda=1.5,kappa=0.2) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  if (is.null(Dshore)){
    Dshore=1/cov_data$ExpShore
  }
  coeff=exp(-kappa*(Dshore/D0)^2)
  omega=omega0/2*(tanh(lambda*(theta+pi/2-atanh(-0.9)/lambda))+tanh(lambda*(theta-(pi/2-atanh(-0.9)/lambda))))*coeff
  return(omega)
}





# Approximation of smooth omega with tensor splines -----------------------

n <- 100000                           #number of observations
D_low=1
D_up=5

theta <- runif(n,-pi,pi)            #sample theta
DistanceShore <- runif(n,D_low,D_up)     #sample DistanceShore

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


#fit with bivariate splines te fo DistanceShore
knots_DistanceShore=list(theta=seq(-pi,pi,len=SP_DF[1]),DistanceShore=seq(D_low,D_up,len=SP_DF[2]))
m1 <- gam(y~te(theta,DistanceShore,k=SP_DF,bs="cs"),knots=knots_DistanceShore)

#visualize
vis.gam(m1);title("tensor product")

df=data.frame("x"=theta,"z"=DistanceShore,"y"=y)
fig <- plot_ly(df,type="scatter3d", x = ~x, y = ~z, z = ~y,mode="markers",marker=list(size=4))
fig



ExpShore=1/DistanceShore

#fit with bivariate splines te of ExpShore
knots_ExpShore=list(theta=seq(-pi,pi,len=SP_DF[1]),ExpShore=seq(1/D_up,1/D_low,len=SP_DF[2]))
m2 <- gam(y~te(theta,ExpShore,k=SP_DF,bs="cs"),knots=knots_ExpShore)


#visualize
vis.gam(m2);title("tensor product")


#get splines coefficients

sp_coeff_Dshore=m1$coefficients
sp_coeff_ExpShore=m2$coefficients


#define spline smooth parameter omega
fomega_splines=function(cov_data) {
  
  if (is.null(cov_data$ExpShore)){
    omega=predict(m1,newdata=cov_data)
  }
  else if (is.null(cov_data$DistanceShore)){
    omega=predict(m2,newdata=cov_data)
  }
  
  return(as.numeric(omega))
}

#k-th basis spline function for omega
fomega_spline_basis=function(cov_data,k=5) {
  
  if (is.null(cov_data$ExpShore)){
    Xp=predict(m1, newdata = cov_data, type = "lpmatrix")
    
  }
  else if (is.null(cov_data$DistanceShore)){
    Xp=predict(m2, newdata = cov_data, type = "lpmatrix")
  }
  coeffs=rep(0,ncol(Xp))
  coeffs[k]=1
  return(as.numeric(Xp%*%coeffs))
}

