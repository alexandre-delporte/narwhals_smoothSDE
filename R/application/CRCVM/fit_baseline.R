
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




#### PREPROCESS EXPLSHORE ----------

D_low=0.07
D_up=3
dataBE1[dataBE1$DistanceShore> D_up,"ExpShore"]=0
dataBE1[dataBE1$DistanceShore< D_low,"ExpShore"]=1/D_low

dataBE2[dataBE2$DistanceShore> D_up,"ExpShore"]=0
dataBE2[dataBE2$DistanceShore< D_low,"ExpShore"]=1/D_low

dataBE1$indicator <- ifelse(dataBE1$DistanceShore < 1, 1, 0)
dataBE2$indicator <- ifelse(dataBE2$DistanceShore < 1, 1, 0)

#### CONSTANT PARAMETERS  -----------------------------------------  

#initial parameters
par0 <- c(0,0,1,4,0)

# Measurement error
sigma_obs=0.035
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))



#model formula
formulas <- list(mu1=~1,mu2=~1,
                 tau =~1,
                 nu=~1,
                 omega=~1)

baseline0<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",response = c("x","y"),
                    par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))

#fit_model
baseline0$fit()
estimates_bas0=as.list(baseline0$tmb_rep(),what="Est")
std_bas0=as.list(baseline0$tmb_rep(),what="Std")

#### ANGLE NORMAL IN OMEGA  -----------------------------------------  

#initial parameters
par0 <- c(0,0,1,4.5,0)

# Measurement error
sigma_obs=0.035
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))



#model formula
formulas <- list(mu1=~1,mu2=~1,
                 tau =~s(ID,bs="re"),
                 nu=~s(ID,bs="re"),
                 omega=~s(AngleNormal,k=4,bs="cs"))

baseline1<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",response = c("x","y"),
                    par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))

#fit_model
baseline1$fit()
estimates_bas1=as.list(baseline1$tmb_rep(),what="Est")
std_bas1=as.list(baseline1$tmb_rep(),what="Std")

res=baseline1$get_all_plots(baseline=NULL,model_name="baseline1",show_CI="pointwise",save=TRUE)

#we see a slight effect of the angle on omega


###  DISTANCE SHORE AND ANGLE NORMAL IN OMEGA  -----------------------

#initial parameters
par0 <- c(0,0,1,4,0)

# Measurement error
sigma_obs=0.05
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

formulas <- list(mu1 = ~1 ,mu2 =~1,tau =~s(ID,bs="re"),nu=~s(ID,bs="re"),
                 omega=~ti(DistanceShore,k=5,bs="cs")+ti(AngleNormal,k=5,bs="cs")+ti(DistanceShore,AngleNormal,k=c(5,5),bs="cs")+s(ID,bs="re"))

baseline2<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),
                    fixpar=c("mu1","mu2"))

#initialize coefficients based on estimates of model3 with all the data
init_coeff_re=model2$coeff_re()
baseline2$update_coeff_re(init_coeff_re)

#fit_model
baseline2$fit()
estimates_bas2=as.list(baseline2$tmb_rep(),what="Est")
std_bas2=as.list(baseline2$tmb_rep(),what="Std")


###  EXP SHORE AND ANGLE NORMAL IN OMEGA  -----------------------

#initial parameters
par0 <- c(0,0,1,4,0)

# Measurement error
sigma_obs=0.05
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

formulas <- list(mu1 = ~1 ,mu2 =~1,tau =~s(ID,bs="re"),nu=~s(ID,bs="re"),
                   omega=~ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,AngleNormal,k=c(5,5),bs="cs")+s(ID,bs="re"))

baseline3<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",
                    response = c("x","y"),par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),
                    fixpar=c("mu1","mu2"))

#initialize coefficients based on estimates of model3 with all the data
init_coeff_re=model3$coeff_re()
baseline3$update_coeff_re(init_coeff_re)

#fit_model
baseline3$fit()
estimates_bas3=as.list(baseline3$tmb_rep(),what="Est")
std_bas3=as.list(baseline3$tmb_rep(),what="Std")



#plot parameters

xmin=list("ExpShore"=1/D_up)
xmax=list("ExpShore"=1/D_low)
link=list("ExpShore"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore")

res=baseline3$get_all_plots(baseline=NULL,model_name="baseline3",
                            xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


##########  ANGLE NORMAL IN OMEGA WITH DISTANCESHORE THRESHOLD


#initial parameters
par0 <- c(0,0,1,4.5,0)

# Measurement error
sigma_obs=0.035
n_obs=length(dataBE2$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))



#model formula
formulas <- list(mu1=~1,mu2=~1,
                 tau =~s(ID,bs="re"),
                 nu=~s(ID,bs="re"),
                 omega=~s(AngleNormal,by=indicator,k=4,bs="cs"))

baseline4<- SDE$new(formulas = formulas,data = dataBE2,type = "RACVM",response = c("x","y"),
                    par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))

#fit_model
baseline4$fit()
estimates_bas4=as.list(baseline4$tmb_rep(),what="Est")
std_bas4=as.list(baseline4$tmb_rep(),what="Std")

res=baseline4$get_all_plots(baseline=NULL,model_name="baseline4",show_CI="pointwise",save=TRUE)


###  AICs VALUES FOR THE BASELINE MODELS-------------------------------

AIC_bas1=baseline1$AIC_marginal()
AIC_bas2=baseline2$AIC_marginal()

