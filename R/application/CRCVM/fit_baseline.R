
# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-06
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/application/constraints/fit_baseline.R
#
# Script Description: fit baseline CRCVM on data before exposure
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space
#SEED FOR REPRODUCTIBILITY
set.seed(42)

library(smoothSDE)
library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(mgcv)



### OUtput txt file --------
filename=paste("fit_baseline_output.txt",sep="")

#Delete file if it exists
if (file.exists(filename)) {
  file.remove(filename)
}

file.create(filename)

##### Get narwhal data ---------

# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
narwhal_data_path <- file.path(par_dir,"Data", "Narwhals")  

# DATA BEFORE EXPOSURE
dataBE1=read.csv(file.path(narwhal_data_path,"DataBE1.csv"), header = TRUE,dec = ".")


cat("Extracting trajectories before exposure and 1 day after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE1[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)


#  DATA BEFORE EXPOSURE WITH 12H  TAGGING EFFECT

dataBE2=read.csv(file.path(narwhal_data_path,"DataBE2.csv"), header = TRUE,dec = ".")

cat("Extracting all trajectories before exposure and 12 hours after tagging ... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE2[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)




#### Preprocess data ----------

#dataBE1=dataBE1[dataBE1$DistanceShore> 0.02,]
#dataBE2=dataBE2[dataBE2$DistanceShore> 0.02,]

#dataBE1[dataBE1$DistanceShore> 3,"ExpShore"]=0
#dataBE1[dataBE1$DistanceShore> 3,"ExpShore"]=0

#### Set measurement errors -----

sigma_obs=0.03
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

### Spline degree of freedom for omega
 SP_DF=7

### Smooth function for omega

fomega=function(cov_data,D0=0.1,omega0=60*pi,lambda=2,kappa=0.2) {
  Dshore=cov_data$DistanceShore
  theta=cov_data$theta
  if (is.null(Dshore)){
    Dshore=1/cov_data$ExpShore
  }
  coeff=exp(-kappa*(Dshore/D0)^2)
  omega=omega0/2*(tanh(lambda*(theta+pi/2-atanh(-0.9)/lambda))+tanh(lambda*(theta-(pi/2-atanh(-0.9)/lambda))))*coeff
  return(omega)
}

#number of points
n <- 1000

#sample theta and DistanceShore points in the domain
theta <- runif(n,-pi,pi);
DistanceShore <- runif(n,0.2,3);
samples=data.frame(theta=theta,DistanceShore=DistanceShore)

#define grid of values
theta_v <- seq(-pi,pi,length=30)
Dshore_v<- seq(0.2,3,length=30)
pr <- data.frame(theta=rep(theta_v,30),DistanceShore=rep(Dshore_v,rep(30,30)))

# true values of the function over this grid
truth <- matrix(fomega(pr),30,30)

#points on the surface perturbed by gaussian noise
f <- fomega(samples)
y <- f+0*rnorm(n)

#plot true function
persp(theta_v,Dshore_v,truth);
title("truth")


#fit with bivariate splines te
m1 <- gam(y~te(theta,DistanceShore,k=SP_DF))

ExpShore=1/DistanceShore

#fit with bivariate splines te
m2 <- gam(y~te(theta,ExpShore,k=SP_DF))

sp_coeff_Dshore=m1$coefficients
sp_coeff_ExpShore=m2$coefficients

names(sp_coeff_Dshore)=paste("omega.",names(sp_coeff_Dshore),sep="")
names(sp_coeff_ExpShore)=paste("omega.",names(sp_coeff_ExpShore),sep="")

fomega_splines=function(cov_data) {
  
  if (is.null(cov_data$ExpShore)){
    omega=predict(m1,newdata=cov_data)
  }
  else if (is.null(cov_data$DistanceShore)){
    omega=predict(m2,newdata=cov_data)
  }
  
  return(as.numeric(omega))
}





#########################  ANGLE NORMAL IN OMEGA   ##############################

#initial parameters
par0 <- c(0,0,1,4,0)

#model formula
formulas <- list(mu1=~1,mu2=~1,
                 tau =~s(AngleNormal,k=6,bs="cs")+s(ID,bs="re"),
                 nu=~s(ID,bs="re"),
                 omega=~s(AngleNormal,k=4,bs="cs"))

baseline1<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",response = c("x","y"),
                    par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))

#fit_model
baseline1$fit()
estimates1=as.list(baseline1$tmb_rep(),what="Est")
std1=as.list(baseline1$tmb_rep(),what="Std")

res=baseline1$get_all_plots(baseline=NULL,model_name="baseline1",show_CI="pointwise",save=TRUE)



#########################  TE SPLINES ANGLE NORMAL AND EXP SHORE IN OMEGA   ##############################

formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~s(AngleNormal,k=6,bs="cs")+s(DistanceShore,k=4,bs="cs"),
                 nu=~1,omega=~te(AngleNormal,DistanceShore,k=c(5,4),bs="cs"))

baseline2<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),
                    fixpar=c("mu1","mu2"))


#fit_model
baseline2$fit()
estimates2=as.list(baseline2$tmb_rep(),what="Est")
std2=as.list(baseline2$tmb_rep(),what="Std")


#plot parameters

xmin=list("DistanceShore"=0.05,"AngleNormal"=-pi)
xmax=list("DistanceShore"=2,"AngleNormal"=pi)
link=list("ExpShore"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore")

res=baseline2$get_all_plots(baseline=NULL,model_name="baseline2",
                            xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)

#########################  TE SPLINES ANGLE NORMAL AND EXP SHORE IN OMEGA WITH FIXED COEFFS   ##############################

formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~s(AngleNormal,k=6,bs="cs")+s(ID,bs="re"),
                 nu=~s(ID,bs="re"),omega=~te(AngleNormal,DistanceShore,k=c(4,3),bs="cs"))

baseline2<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",
                    response = c("x","y"),par0 = par0,other_data=list("H"=H),
                    fixpar=c("mu1","mu2"),map=list(coeff_re=factor(c(1:6,rep(NA,SP_DF^2-1))),coeff_fe=factor(c(NA,NA,7:8,NA))))

baseline2$update_coeff_re(c(rep(0,6),sp_coeff_Dshore[-1]))
new_coeff_fe=baseline2$coeff_fe()
new_coeff_fe[5]=sp_coeff_Dshore[1]
baseline2$update_coeff_fe(new_coeff_fe)

#fit_model
baseline2$fit()
estimates2=as.list(baseline2$tmb_rep(),what="Est")
std2=as.list(baseline2$tmb_rep(),what="Std")


#plot parameters

xmin=list("DistanceShore"=0.05,"AngleNormal"=-pi)
xmax=list("DistanceShore"=2,"AngleNormal"=pi)
link=list("ExpShore"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore")

res=baseline2$get_all_plots(baseline=NULL,model_name="baseline2",
                            xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)

#########################  TI SPLINES ANGLE NORMAL AND EXPSHORE IN OMEGA  ##############################

formulas <- list(mu1 = ~1 ,mu2 =~1,tau = ~1,
                 nu=~1,omega=~te(AngleNormal,ExpShore,k=3,bs="cs"))

baseline3<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(0.04)),
                    fixpar=c("mu1","mu2"))

#fit_model
baseline3$fit()
estimates3=as.list(baseline3$tmb_rep(),what="Est")
std3=as.list(baseline3$tmb_rep(),what="Std")


res=baseline3$get_all_plots(baseline=NULL,model_name="baseline3",
                            xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


###########################         AICs VALUES FOR THE BASELINE MODELS        ####################################

AIC_bas1=baseline1$AIC_marginal()
AIC_bas2=baseline2$AIC_marginal()
AIC_bas3=baseline3$AIC_marginal()



cat("--------------------- BASELINE MODELS ---------------------: \n",
    baseline1$formulas(),AIC_bas1,
    baseline2$formulas(),AIC_bas2,baseline3$formulas(),AIC_bas3,file=filename,sep="\n",append=TRUE)
