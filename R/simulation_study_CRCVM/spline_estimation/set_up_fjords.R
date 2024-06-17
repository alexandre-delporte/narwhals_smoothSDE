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
library(smoothSDE)          #to compute nearest shore point

#https://stats.stackexchange.com/questions/37647/what-is-the-minimum-recommended-number-of-groups-for-a-random-effects-factor
#It is recommended to have around 10 groups to get good estimates of the random effect variance
N_ID=6                   #number of individual tracks per batch
TMAX=24*5                  #duration of each track in hour
SP_DF=c(3,3)                #degree of freedom in tensor splines
TAU_0=1.5                      # tau intercept
NU_0=4                       #nu intercept
SIGMA_TAU=0.2                #variance of random effects on tau
SIGMA_NU=0.1                 # variance of random effects on nu
DMIN=1                      #min distance to shore for initial position
DELTA=1/60                   #time step
PAR0=c(0,0,1,1,0)            # initial values (mu1,mu2,tau,nu,omega)
SIGMA_OBS=0.05             # measurement error
BY=10


# Land polygons --------------------


# Set the path to the directory containing the data
dir=dirname(dirname(dirname(dirname(getwd())))) #directory of project
greenland_data_path <- file.path(dir,"Data", "Greenland")

#get the land and coastline geometries from the geojson file
border<-st_read(file.path(greenland_data_path,"land_utm.shp"))
st_crs(border)="+init=EPSG:32626"
border<-st_make_valid(border)

v0=c(0,0)

#generate random initial points
x0=matrix(rep(NA,N_ID*2),ncol=2)
colnames(x0)=c("x1","x2")

i=1
while (i<=N_ID) {
  #choose location uniformly in the map
  x=c(runif(1,min=430,max=500),runif(1,min=7760,max=7900))
  p=nearest_shore_point(st_point(x),border)
  Dshore=(x[1]-p[1])^2+(x[2]-p[2])^2
  #keep it as initial position if it is at least 50 metres away from the shore
  if (Dshore>DMIN) {
    x0[i,]=x
    i=i+1
  }
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

fomega=function(cov_data,D0=0.5,omega0=60*pi,lambda=1.5,kappa=0.2) {
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


#defune spline smooth parameter omega
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
fomega_spline_basis=function(cov_data,k=1) {
  
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

