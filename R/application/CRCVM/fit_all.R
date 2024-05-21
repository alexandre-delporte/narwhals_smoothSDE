
# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-16
#
# Script Name:    ~/Documents/Research/Projects/narwhals_smoothSDE/R/application/CRCVM/fit_all.R
#
# Script Description: fit sde models with shore effects on all the data (before + after exposure).
# The spline coefficients estimated here are used as initial values for the baseline model.
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console


set.seed(42)

library(smoothSDE)
library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)


#   GET NARWHAL DATA   ---------------------------------


# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
narwhal_data_path <- file.path(par_dir,"Data", "Narwhals")  


# DATA BEFORE EXPOSURE

allData=read.csv(file.path(narwhal_data_path,"allData.csv"), header = TRUE,dec = ".")




allData=allData[allData$time>12,]

###SOME PREPROCESSING FOR EXPSHORE ------------------------------------

#treshold on distance to shore
D_low=0.07
D_up=3

allData[allData$DistanceShore> D_up,"ExpShore"]=0
allData[allData$DistanceShore< D_low,"ExpShore"]=1/D_low




# MODEL WITH ANGLENORMAL IN OMEGA  -----------------------------

#initial parameters
par0 <- c(0,0,2,4.5,0)

#measurement error
sigma_obs=0.05
n_obs=length(allData$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))


#model formula
formulas <- list(mu1=~1,mu2=~1,tau =~s(ID,bs="re"),nu=~s(ID,bs="re"),
                 omega=~s(AngleNormal,k=4,bs="cs"))

model1<- SDE$new(formulas = formulas,data = allData[allData$time<2.5*24,],type = "RACVM",response = c("x","y"),
                    par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),fixpar=c("mu1","mu2"))

#fit_model
model1$fit()
estimates_mod1=as.list(model1$tmb_rep(),what="Est")
std_mod1=as.list(model1$tmb_rep(),what="Std")

res=model1$get_all_plots(baseline=NULL,model_name="model1",show_CI="pointwise",save=TRUE)

#we see a slight effect of the angle on omega 

##   MODEL WITH ANGLENORMAL AND DISTANCE SHORE IN OMEGA -------------------

#initial parameters
par0 <- c(0,0,1,4,0)

#measurement error
sigma_obs=0.05

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(ID,bs="re"),nu=~s(ID,bs="re"),
                 omega=~ti(DistanceShore,k=5,bs="cs")+ti(AngleNormal,k=5,bs="cs")+ti(DistanceShore,AngleNormal,k=5,bs="cs")+s(ID,bs="re"))

model2<- SDE$new(formulas = formulas,data = allData,type = "RACVM",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),
                    other_data=list("log_sigma_obs0"=log(sigma_obs)))



#fit_model
model2$fit()
estimates_mod2=as.list(model2$tmb_rep(),what="Est")
std_mod2=as.list(model1$tmb_rep(),what="Std")


#plot parameters


xmin=list("DistanceShore"=D_low,"AngleNormal"=-pi)
xmax=list("DistanceShore"=D_up,"AngleNormal"=pi)
xlabel=list("DistanceShore"="Distance to shore")
res=model2$get_all_plots(baseline=NULL,model_name="model2",xmin=xmin,xmax=xmax,xlabel=xlabel,show_CI="pointwise",save=TRUE)



#########################  RACVM  MODEL WITH tensor splines of  AngleNormal, ExpShore in omega ##############################

#initial parameters
par0 <- c(0,0,1,4,0)

#measurement error
sigma_obs=0.05

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(ID,bs="re"),nu=~s(ID,bs="re"),
                 omega=~ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,AngleNormal,k=c(5,5),bs="cs")+s(ID,bs="re"))

model3<- SDE$new(formulas = formulas,data =allData,type = "RACVM",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),
                 other_data=list("log_sigma_obs0"=log(sigma_obs)))



#fit_model
model3$fit()
estimates_mod3=as.list(model3$tmb_rep(),what="Est")
std_mod3=as.list(model3$tmb_rep(),what="Std")

xmin=list("ExpShore"=1/D_up,"AngleNormal"=-pi)
xmax=list("ExpShore"=1/D_low,"AngleNormal"=pi)
link=list("ExpShore"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore")

#plot parameters

model3$get_all_plots(baseline=NULL,model_name="model3",xmin=xmin,
                        xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)



#########################   Model with AngleNormal, ExpShore in omega and AngleNormal in tau ##############################

#initial parameters
par0 <- c(0,0,1,4,0)

#measurement error
sigma_obs=0.05

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(AngleNormal,k=5,bs="cs")+s(ID,bs="re"),
                 nu=~s(ID,bs="re"),
                 omega=~ti(AngleNormal,k=4,bs="cs")+ti(ExpShore,k=4,bs="cs")+ti(AngleNormal,ExpShore,k=4,bs="cs")+s(ID,bs="re"))

model4<- SDE$new(formulas = formulas,data = allData,type = "RACVM",
                 response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),other_data=list("log_sigma_obs0"=log(sigma_obs)))



#fit_model
model4$fit()
estimates_mod4=as.list(model4$tmb_rep(),what="Est")
std_mod4=as.list(model4$tmb_rep(),what="Std")

#plot parameters
xmin=list("ExpShore"=1/D_up,"AngleNormal"=-pi,"ExpShip"=1/45)
xmax=list("ExpShore"=1/D_low,"AngleNormal"=pi,"ExpShip"=1/3)
link=list("ExpShore"=(\(x) 1/x),"ExpShip"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore","ExpShip"="Distance to ship")

model4$get_all_plots(baseline=NULL,model_name="model4",xmin=xmin,
                     xmax=xmax,link=link,xlabel=xlabel,show_CI="none",save=TRUE)



#########################  Model with AngleNormal, ExpShore in omega and ExpShip in nu ##############################

#initial parameters
par0 <- c(0,0,1,4,0)

#measurement error
sigma_obs=0.05

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(ID,bs="re"),nu=~s(ExpShip,k=3,bs="cs")+s(ID,bs="re"),
                 omega=~ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,ExpShore,k=5,bs="cs")+s(ID,bs="re"))

model5<- SDE$new(formulas = formulas,data = allData,type = "RACVM",response = c("x","y"),par0 = par0,
                 fixpar=c("mu1","mu2"),other_data=list("log_sigma_obs0"=log(sigma_obs)))


#fit_model
model5$fit()
estimates_mod5=as.list(model5$tmb_rep(),what="Est")
std_mod5=as.list(model5$tmb_rep(),what="Std")

#plot parameters
xmin=list("ExpShore"=1/D_up,"AngleNormal"=-pi,"ExpShip"=1/80)
xmax=list("ExpShore"=1/D_low,"AngleNormal"=pi,"ExpShip"=1/3)
link=list("ExpShore"=(\(x) 1/x),"ExpShip"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore","ExpShip"="Distance to ship")


model5$get_all_plots(baseline=NULL,model_name="model5",xmin=xmin,
              xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)



#########################   Model with AngleNormal, ExpShore in omega and ExpShip in nu and tau ##############################

#initial parameters
par0 <- c(0,0,1,4,0)

#measurement error
sigma_obs=0.05

#define model
formulas <- list(mu1=~1,mu2=~1,tau =~s(ExpShip,k=3,bs="cs")+s(ID,bs="re"),
                 nu=~s(ExpShip,k=3,bs="cs")+s(ID,bs="re"),
                 omega=~ti(AngleNormal,k=5,bs="cs")+ti(ExpShore,k=5,bs="cs")+ti(AngleNormal,ExpShore,k=5,bs="cs")+s(ID,bs="re"))

model6<- SDE$new(formulas = formulas,data = allData,type = "RACVM",
                    response = c("x","y"),par0 = par0,fixpar=c("mu1","mu2"),other_data=list("log_sigma_obs0"=log(sigma_obs)))



#fit_model
model6$fit()
estimates_mod6=as.list(model6$tmb_rep(),what="Est")
std_mod=as.list(model6$tmb_rep(),what="Std")

#plot parameters
xmin=list("ExpShore"=1/D_up,"AngleNormal"=-pi,"ExpShip"=1/45)
xmax=list("ExpShore"=1/D_low,"AngleNormal"=pi,"ExpShip"=1/3)
link=list("ExpShore"=(\(x) 1/x),"ExpShip"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore","ExpShip"="Distance to ship")

model6$get_all_plots(baseline=NULL,model_name="model6",xmin=xmin,
                        xmax=xmax,link=link,xlabel=xlabel,show_CI="none",save=TRUE)




#############################      model AIC VALUES    #####################################################

AIC1=model1$AIC_marginal()
AIC2=model2$AIC_marginal()
AIC3=model3$AIC_marginal()
AIC4=model4$AIC_marginal()
AIC5=model5$AIC_marginal()


